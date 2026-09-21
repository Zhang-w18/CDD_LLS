from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PDCCHDecodeResult:
    payload_bits: np.ndarray
    crc_status: np.ndarray


def rnti_crc_mask(rnti: int) -> np.ndarray:
    value = int(rnti)
    if value < 0 or value > 0xFFFF:
        raise ValueError("RNTI must lie in [0,0xFFFF].")
    bits = np.asarray([(value >> (15 - i)) & 1 for i in range(16)], dtype=np.int8)
    return np.concatenate([np.zeros(8, dtype=np.int8), bits])


def dci_mother_code_length(payload_bits: int, coded_bits: int) -> int:
    """Return 38.212 DCI Polar mother-code length N."""
    information_bits = int(payload_bits) + 24
    target_bits = int(coded_bits)
    ceiling_log2 = int(np.ceil(np.log2(target_bits)))
    if (
        target_bits <= (9.0 / 8.0) * 2 ** (ceiling_log2 - 1)
        and information_bits / target_bits < 9.0 / 16.0
    ):
        n1 = ceiling_log2 - 1
    else:
        n1 = ceiling_log2
    n2 = int(np.ceil(np.log2(8 * information_bits)))
    return int(2 ** max(min(n1, n2, 9), 5))


class SionnaPDCCHPolarCodec:
    """3GPP DCI wrapper around Sionna's downlink Polar implementation.

    Sionna 1.0.2 intentionally omits the DCI all-one CRC initialization and
    RNTI parity-bit mask. Both are deterministic affine offsets. This wrapper
    propagates the exact offset through input interleaving, Polar encoding and
    rate matching, applies it at the transmitter, and removes it in the LLR
    domain before Sionna's CRC-aided decoder.
    """

    def __init__(
        self,
        payload_bits: int,
        coded_bits: int,
        crc_rnti: int,
        scrambling_id: int = 0,
        data_scrambling_rnti: int = 0,
        list_size: int = 8,
        decoder_type: str = "SCL",
        cpu_only: bool = True,
    ) -> None:
        import tensorflow as tf
        from sionna.phy.fec.polar import (
            Polar5GDecoder,
            Polar5GEncoder,
            PolarEncoder,
            PolarSCLDecoder,
        )
        from sionna.phy.nr import generate_prng_seq

        self.tf = tf
        self.payload_bits = int(payload_bits)
        self.coded_bits = int(coded_bits)
        self.crc_rnti = int(crc_rnti)
        self.scrambling_id = int(scrambling_id)
        self.data_scrambling_rnti = int(data_scrambling_rnti)
        self.list_size = int(list_size)
        self.decoder_type = str(decoder_type)

        if self.payload_bits <= 0 or self.payload_bits > 140:
            raise ValueError("DCI payload_bits must lie in [1,140].")
        if self.coded_bits != int(self.coded_bits) or self.coded_bits <= self.payload_bits + 24:
            raise ValueError("coded_bits must exceed payload_bits+24.")
        if self.list_size <= 0:
            raise ValueError("list_size must be positive.")
        rnti_crc_mask(self.crc_rnti)
        if self.scrambling_id < 0 or self.scrambling_id > 65535:
            raise ValueError("scrambling_id must lie in [0,65535].")
        if self.data_scrambling_rnti < 0 or self.data_scrambling_rnti > 0xFFFF:
            raise ValueError("data_scrambling_rnti must lie in [0,0xFFFF].")

        self.mother_code_length = dci_mother_code_length(
            self.payload_bits, self.coded_bits
        )
        self.external_repetition = self.coded_bits > 576
        if self.external_repetition and self.coded_bits <= self.mother_code_length:
            raise ValueError("External rate matching requires coded_bits>N.")
        self.sionna_target_bits = (
            self.mother_code_length if self.external_repetition else self.coded_bits
        )
        if self.external_repetition:
            self.rate_matching_mode = "repetition"
        elif self.coded_bits < self.mother_code_length:
            self.rate_matching_mode = (
                "puncturing"
                if (self.payload_bits + 24) / self.coded_bits <= 7.0 / 16.0
                else "shortening"
            )
        else:
            self.rate_matching_mode = "repetition" if self.coded_bits > self.mother_code_length else "none"

        self.encoder = Polar5GEncoder(
            self.payload_bits,
            self.sionna_target_bits,
            channel_type="downlink",
            precision="double",
        )
        self.decoder = Polar5GDecoder(
            self.encoder,
            dec_type=self.decoder_type,
            list_size=self.list_size,
            return_crc_status=True,
            precision="double",
        )
        if bool(cpu_only) and self.decoder_type in ("SCL", "hybSCL"):
            # Polar5GDecoder 1.0.2 does not forward cpu_only to its SCL layer.
            # Replacing that layer keeps the public 5G rate-recovery wrapper
            # while making the documented CPU option effective.
            self.decoder._polar_dec = PolarSCLDecoder(
                self.encoder.frozen_pos,
                self.encoder.n_polar,
                crc_degree=self.encoder.enc_crc.crc_degree,
                list_size=self.list_size,
                use_hybrid_sc=self.decoder_type == "hybSCL",
                cpu_only=True,
                ind_iil_inv=self.decoder.ind_iil_inv,
                precision="double",
            )
            self.decoder._dec_crc = self.decoder._polar_dec._crc_decoder

        parity_offset = self._build_crc_parity_offset()
        information_offset = np.concatenate(
            [np.zeros(self.payload_bits, dtype=np.int8), parity_offset]
        )
        interleaved = information_offset[np.asarray(self.encoder._ind_input_int, dtype=np.int64)]
        offset_encoder = PolarEncoder(
            self.encoder.frozen_pos,
            self.encoder.n_polar,
            precision="double",
        )
        mother = offset_encoder(
            tf.constant(interleaved[None, :], dtype=tf.float64)
        ).numpy().astype(np.int8)[0]
        base_codeword_offset = mother[
            np.asarray(self.encoder._ind_rate_matching, dtype=np.int64)
        ].astype(np.int8)
        self.codeword_offset = self._expand_external_repetition(base_codeword_offset)

        c_init = (
            (self.data_scrambling_rnti << 16) + self.scrambling_id
        ) % (1 << 31)
        self.scrambling_sequence = generate_prng_seq(
            self.coded_bits, int(c_init)
        ).astype(np.int8)

    def _expand_external_repetition(self, values: np.ndarray) -> np.ndarray:
        arr = np.asarray(values)
        if not self.external_repetition:
            return arr
        indices = np.arange(self.coded_bits, dtype=np.int64) % self.sionna_target_bits
        return arr[..., indices]

    def _recover_external_repetition(self, values: np.ndarray) -> np.ndarray:
        arr = np.asarray(values)
        if not self.external_repetition:
            return arr
        quotient, remainder = divmod(self.coded_bits, self.sionna_target_bits)
        recovered = arr[..., : quotient * self.sionna_target_bits].reshape(
            *arr.shape[:-1], quotient, self.sionna_target_bits
        ).sum(axis=-2)
        if remainder:
            recovered[..., :remainder] += arr[..., quotient * self.sionna_target_bits :]
        return recovered

    def _crc_encode(self, bits: np.ndarray) -> np.ndarray:
        arr = np.asarray(bits, dtype=np.int8)
        if arr.ndim == 1:
            arr = arr[None, :]
        encoded = self.encoder.enc_crc(
            self.tf.constant(arr, dtype=self.tf.float64)
        ).numpy().astype(np.int8)
        return encoded

    def _build_crc_parity_offset(self) -> np.ndarray:
        zeros = np.zeros(self.payload_bits, dtype=np.int8)
        base_crc = self._crc_encode(zeros)[0, -24:]
        prefixed = np.concatenate([np.ones(24, dtype=np.int8), zeros])
        initialized_crc = self._crc_encode(prefixed)[0, -24:]
        return (
            base_crc ^ initialized_crc ^ rnti_crc_mask(self.crc_rnti)
        ).astype(np.int8)

    def standard_crc_bits(self, payload: np.ndarray) -> np.ndarray:
        arr = np.asarray(payload, dtype=np.int8)
        if arr.ndim == 1:
            arr = arr[None, :]
        if arr.shape[-1] != self.payload_bits:
            raise ValueError("Payload width does not match codec payload_bits.")
        prefix = np.ones((arr.shape[0], 24), dtype=np.int8)
        initialized = self._crc_encode(np.concatenate([prefix, arr], axis=1))[:, -24:]
        return (initialized ^ rnti_crc_mask(self.crc_rnti)[None, :]).astype(np.int8)

    def encode_unscrambled(self, payload: np.ndarray) -> np.ndarray:
        arr = np.asarray(payload, dtype=np.int8)
        if arr.ndim == 1:
            arr = arr[None, :]
        if arr.shape[-1] != self.payload_bits:
            raise ValueError("Payload width does not match codec payload_bits.")
        encoded = self.encoder(
            self.tf.constant(arr, dtype=self.tf.float64)
        ).numpy().astype(np.int8)
        encoded = self._expand_external_repetition(encoded)
        return (encoded ^ self.codeword_offset[None, :]).astype(np.int8)

    def encode(self, payload: np.ndarray) -> np.ndarray:
        exact = self.encode_unscrambled(payload)
        return (exact ^ self.scrambling_sequence[None, :]).astype(np.int8)

    def prepare_decoder_llr(self, received_llr: np.ndarray) -> np.ndarray:
        llr = np.asarray(received_llr, dtype=np.float64)
        if llr.shape[-1] != self.coded_bits:
            raise ValueError("LLR width does not match codec coded_bits.")
        flips = (
            self.scrambling_sequence ^ self.codeword_offset
        ).astype(np.int8)
        prepared = llr * (1.0 - 2.0 * flips[None, :])
        return self._recover_external_repetition(prepared)

    def decode(self, received_llr: np.ndarray) -> PDCCHDecodeResult:
        prepared = self.prepare_decoder_llr(received_llr)
        payload, crc_status = self.decoder(
            self.tf.constant(prepared, dtype=self.tf.float64)
        )
        return PDCCHDecodeResult(
            payload_bits=np.rint(payload.numpy()).astype(np.int8),
            crc_status=np.asarray(crc_status.numpy(), dtype=bool).reshape(-1),
        )

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence
import os
import numpy as np


@dataclass(frozen=True)
class LDPCDecodeResult:
    decoded_bits_by_cb: List[np.ndarray]
    cb_success: List[bool]
    tb_success: bool
    cb_bler: float
    tb_bler: float
    goodput_bits: int


class SionnaLDPCAdapter:
    """Small Sionna LDPC5G adapter for one transport block split into CBs."""

    def __init__(
        self,
        cb_k_values: Sequence[int],
        cb_e_values: Sequence[int],
        num_iter: int = 20,
        llr_clip: float = 50.0,
    ):
        if len(cb_k_values) != len(cb_e_values):
            raise ValueError("cb_k_values and cb_e_values must have the same length.")
        self.cb_k_values = [int(x) for x in cb_k_values]
        self.cb_e_values = [int(x) for x in cb_e_values]
        self.num_iter = int(num_iter)
        self.llr_clip = float(llr_clip)
        self.tf, self.LDPC5GEncoder, self.LDPC5GDecoder = self._import_ldpc()
        self.encoder_by_cb = {}
        self.decoder_by_cb = {}
        for cb_idx, (k, n) in enumerate(zip(self.cb_k_values, self.cb_e_values)):
            enc = self.LDPC5GEncoder(k=int(k), n=int(n))
            dec = self._make_decoder(enc)
            self.encoder_by_cb[int(cb_idx)] = enc
            self.decoder_by_cb[int(cb_idx)] = dec

    @staticmethod
    def _import_ldpc():
        os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
        os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/cdd_lls_matplotlib")
        os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
        import tensorflow as tf

        try:
            from sionna.phy.fec.ldpc import LDPC5GDecoder, LDPC5GEncoder
        except Exception:
            from sionna.fec.ldpc.decoding import LDPC5GDecoder
            from sionna.fec.ldpc.encoding import LDPC5GEncoder
        return tf, LDPC5GEncoder, LDPC5GDecoder

    def _make_decoder(self, encoder):
        try:
            return self.LDPC5GDecoder(encoder, hard_out=True, num_iter=self.num_iter)
        except TypeError:
            try:
                return self.LDPC5GDecoder(
                    encoder,
                    hard_out=True,
                    num_iter=self.num_iter,
                    return_infobits=True,
                )
            except TypeError:
                return self.LDPC5GDecoder(encoder, hard_out=True)

    def encode_one(self, cb_index: int, payload_bits: np.ndarray) -> np.ndarray:
        cb_index = int(cb_index)
        bits = np.asarray(payload_bits, dtype=np.float32).reshape(1, -1)
        expected = int(self.cb_k_values[cb_index])
        if bits.shape[1] != expected:
            raise ValueError(f"CB{cb_index} payload length {bits.shape[1]} does not match k={expected}.")
        c = self.encoder_by_cb[cb_index](self.tf.constant(bits, dtype=self.tf.float32))
        return np.rint(c.numpy()[0]).astype(np.int8)

    def encode(self, payload_bits_by_cb: Sequence[np.ndarray]) -> List[np.ndarray]:
        if len(payload_bits_by_cb) != len(self.cb_k_values):
            raise ValueError("payload_bits_by_cb length does not match adapter CB count.")
        return [self.encode_one(i, bits) for i, bits in enumerate(payload_bits_by_cb)]

    def decode_one(self, cb_index: int, llr: np.ndarray) -> np.ndarray:
        cb_index = int(cb_index)
        arr = np.asarray(llr, dtype=np.float32).reshape(1, -1)
        expected = int(self.cb_e_values[cb_index])
        if arr.shape[1] != expected:
            raise ValueError(f"CB{cb_index} LLR length {arr.shape[1]} does not match E={expected}.")
        if self.llr_clip > 0:
            arr = np.clip(arr, -self.llr_clip, self.llr_clip)
        b_hat = self.decoder_by_cb[cb_index](self.tf.constant(arr, dtype=self.tf.float32))
        return np.rint(b_hat.numpy()[0]).astype(np.int8)

    def decode(
        self,
        llrs_by_cb: Sequence[np.ndarray],
        reference_payload_bits_by_cb: Sequence[np.ndarray],
    ) -> LDPCDecodeResult:
        if len(llrs_by_cb) != len(self.cb_k_values):
            raise ValueError("llrs_by_cb length does not match adapter CB count.")
        decoded = [self.decode_one(i, llr) for i, llr in enumerate(llrs_by_cb)]
        return summarize_decode(decoded, reference_payload_bits_by_cb)


def summarize_decode(
    decoded_bits_by_cb: Sequence[np.ndarray],
    reference_payload_bits_by_cb: Sequence[np.ndarray],
) -> LDPCDecodeResult:
    decoded = [np.asarray(x, dtype=np.int8).reshape(-1) for x in decoded_bits_by_cb]
    refs = [np.asarray(x, dtype=np.int8).reshape(-1) for x in reference_payload_bits_by_cb]
    if len(decoded) != len(refs):
        raise ValueError("decoded and reference CB counts differ.")
    if not decoded:
        raise ValueError("At least one CB is required.")

    cb_success = []
    goodput_bits = 0
    for got, ref in zip(decoded, refs):
        if got.shape != ref.shape:
            raise ValueError("decoded and reference CB shapes differ.")
        ok = bool(np.array_equal(got, ref))
        cb_success.append(ok)
        if ok:
            goodput_bits += int(ref.size)

    tb_success = bool(all(cb_success))
    cb_errors = sum(1 for ok in cb_success if not ok)
    return LDPCDecodeResult(
        decoded_bits_by_cb=decoded,
        cb_success=cb_success,
        tb_success=tb_success,
        cb_bler=float(cb_errors) / float(len(cb_success)),
        tb_bler=0.0 if tb_success else 1.0,
        goodput_bits=int(goodput_bits),
    )


def ldpc_backend_available() -> tuple[bool, str]:
    try:
        SionnaLDPCAdapter([128], [256], num_iter=1)
        return True, "ok"
    except Exception as exc:
        return False, repr(exc)

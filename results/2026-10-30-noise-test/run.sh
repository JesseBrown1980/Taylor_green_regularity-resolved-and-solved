#!/bin/bash
# SH layer (HBI -> HBP -> SHA -> SH -> HASH, json=0). Checks the seals and the pinned trit inputs,
# builds tgrust with Rust 1.81 (std only), runs clippy -D warnings, and checks that the run
# reproduces RESULT.hbp byte for byte. Never overwrites a record.
set -u
cd "$(dirname "$0")"
sha256sum -c --quiet PREREG.md.sha256 ADDENDUM.md.sha256 INIT.sha256 || { echo "TG_REPLAY|result=INCONCLUSIVE|reason=seal_or_input|json=0"; exit 2; }
export CARGO_TARGET_DIR=$(mktemp -d)
cargo +1.81.0 build -q --release --locked --manifest-path tgrust/Cargo.toml || { echo "TG_REPLAY|result=INCONCLUSIVE|reason=build|json=0"; exit 3; }
cargo +1.81.0 clippy -q --release --locked --manifest-path tgrust/Cargo.toml -- -D warnings || { echo "TG_REPLAY|result=FAIL|reason=clippy|json=0"; exit 4; }
same=0; "$CARGO_TARGET_DIR/release/tgrust" INIT_TGV.trits INIT_NULL.trits GATES.lut | cmp -s - RESULT.hbp && same=1
echo "TG_REPLAY|rustc=$(rustc +1.81.0 --version | tr ' ' '_')|byte_identical=$same|json=0"
rm -rf "$CARGO_TARGET_DIR"
[ $same = 1 ]

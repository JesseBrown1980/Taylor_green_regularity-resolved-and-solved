#!/bin/bash
# SH layer (HBI -> HBP -> SHA -> SH -> HASH, json=0): checks seals and pinned inputs, builds tgwind
# with Rust 1.81 (std only), clippy -D warnings, checks byte identity with RESULT.hbp, and checks
# that RESULT.run1.hbp equals RESULT.hbp without the addendum-2 rows.
set -u
cd "$(dirname "$0")"
sha256sum -c --quiet PREREG.md.sha256 ADDENDUM.md.sha256 ADDENDUM2.md.sha256 INIT.sha256 || { echo "TGW_REPLAY|result=INCONCLUSIVE|reason=seal|json=0"; exit 2; }
export CARGO_TARGET_DIR=$(mktemp -d)
cargo +1.81.0 build -q --release --locked --manifest-path tgwind/Cargo.toml || { echo "TGW_REPLAY|result=INCONCLUSIVE|reason=build|json=0"; exit 3; }
cargo +1.81.0 clippy -q --release --locked --manifest-path tgwind/Cargo.toml -- -D warnings || { echo "TGW_REPLAY|result=FAIL|reason=clippy|json=0"; exit 4; }
same=0; "$CARGO_TARGET_DIR/release/tgwind" INIT_TGV.trits GATES.lut | cmp -s - RESULT.hbp && same=1
r1=0; grep -vE '^TGW(ROW2|SUM2|RESULT2)' RESULT.hbp | cmp -s - RESULT.run1.hbp && r1=1
echo "TGW_REPLAY|rustc=$(rustc +1.81.0 --version | tr ' ' '_')|byte_identical=$same|run1_rows_identical=$r1|json=0"
rm -rf "$CARGO_TARGET_DIR"
[ $same = 1 ] && [ $r1 = 1 ]

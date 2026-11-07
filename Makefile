.PHONY: check fmt clippy test doc coverage audit msrv docs-check docs bench

check: fmt clippy test doc coverage audit msrv docs-check

fmt:
	cargo fmt --all --check

clippy:
	cargo clippy --locked --all-targets --all-features -- -D warnings

test:
	cargo test --locked --all-features
	cargo test --locked

doc:
	RUSTDOCFLAGS="-D warnings" cargo doc --locked --no-deps --all-features

coverage:
	cargo llvm-cov --locked --all-features --fail-under-lines 100

audit:
	cargo audit

# MSRV is imposed by geometry_r (let chains): Rust 1.88.
msrv:
	cargo +1.88.0 check --locked --lib --all-features

docs-check:
	@bash scripts/docs/check_docs.sh
	@python3 scripts/docs/check_links.py .

docs: docs-check
	@echo "Documentation checks passed"

bench:
	cargo bench

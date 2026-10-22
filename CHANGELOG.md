# CHANGELOG
<!-- version list -->

## v0.3.0 (2026-09-23)

### Documentation

- Add README badges ([#11](https://github.com/thedannymarsh/nd_line/pull/11), [`a7d62ba`](https://github.com/thedannymarsh/nd_line/commit/a7d62ba66aa01ce15d95dc0d41e94664ec7d27cb))

- Seed CHANGELOG from v0.2.1 release and polish README notes ([#13](https://github.com/thedannymarsh/nd_line/pull/13), [`a920d18`](https://github.com/thedannymarsh/nd_line/commit/a920d18eeca3de950c3193b36129caa3024708fc))

### Features

- Replace splineify with to_spline that returns a new line ([#5](https://github.com/thedannymarsh/nd_line/pull/5), [`65eebce`](https://github.com/thedannymarsh/nd_line/commit/65eebce5188077d8040ecc729ecf3723f3fc52ab))

### Breaking Changes

- Splineify() is now to_spline() and returns a new nd_line instead of mutating in place. major_on_zero is false so 0.x breaking changes bump minor rather than 1.0.0.

## v0.2.1 (2026-09-14)

### Bug Fixes

- Update return values from splprep ([`08903d5`](https://github.com/thedannymarsh/nd_line/commit/08903d576983a781f8d321c0fea4295f1a5070c9))

- Update splprep return values ([`5c02a36`](https://github.com/thedannymarsh/nd_line/commit/5c02a3615286858f263a473f6cece8071cb79be6))

### Build System

- Add python semantic release, pre-commit ci, and upload to pypi ([`2d69cae`](https://github.com/thedannymarsh/nd_line/commit/2d69cae5d0e5291b7dc434137912b6b9697f5dfc))

# v2.5.31a -- validator false-positive fix (caught by the gate on the first real run)

The v2.5.31 gate did its job on the first real ABFE submission: it refused to
submit. But its verdict was WRONG for two patterns that are actually valid, so it
would have blocked good runs. Both are now fixed.

FALSE POSITIVE 1 -- inline "&wt type='END' /":
  _restraint_block writes the END terminator inline in min/dens/eq/prod. pmemd
  TOLERATES the bare END terminator (it carries no istep/value data to mis-read);
  decharge/vdw legs have used it for years. The v2.5.29 failure was the data-bearing
  TEMP_0 card, not END. Fix: only flag inline &wt cards whose type != END (or that
  carry data fields). The bare inline "&wt type='END' /" is allowed.

FALSE POSITIVE 2 -- empty timask2 under icfe=1:
  For ABFE (absolute decoupling of ONE molecule) timask2 is correctly empty; the
  driver logs timask1=':LIG' timask2=''. Requiring a non-empty timask2 there is
  wrong. Fix: require timask1 always; require timask2 ONLY when scmask2 is
  non-empty (i.e. a genuine two-region/RBFE transformation).

STILL CAUGHT (verified): genuine inline TEMP_0 data card (2.5.29 bug); '*' mask
wildcard (2.5.28); missing -ref on a restrained stage (2.5.27); non-NVT heat;
RBFE with scmask2 set but timask2 empty; missing DISANG file.

ACTION: regenerate with 2.5.31a -- the complex_restraint (and all) legs should now
validate clean and submit.

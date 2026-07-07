# v2.5.0 build final42: GUI NameError on Setup_and_Launch page

## Symptom
Opening the Streamlit "0 - Setup and Launch" page crashed immediately for
EVERY method (MM-GBSA / RBFE / ABFE):

    File ".../gui/pages/0_Setup_and_Launch.py", line 159, in <module>
        params["ligand_resname"] = st.text_input(
    NameError: name 'params' is not defined

## Root cause
The ligand-residue-name widget (section "2 - Inputs") writes into `params`
(params["ligand_resname"] = ...), but `params` was only initialized later in
section "3 - Parameters" (`params: dict = {}`), ~80 lines further down. On
first render the global `params` did not yet exist -> NameError. (The line-78
`params: dict = field(default_factory=dict)` is an unrelated dataclass FIELD,
not the module global.)

## Fix
* Initialize `params: dict = {}` in section 2, immediately after the
  "2 - Inputs" subheader, so early writers (ligand_resname) are safe.
* Section 3 no longer unconditionally re-binds `params` (which would have
  wiped ligand_resname); it now only creates it if missing:
      try: params
      except NameError: params = {}

## Scope
GUI only. No effect on the CLI driver, the running ABFE production job, or the
analysis fixes (final38-41). MM-GBSA / RBFE / ABFE setup pages all load now.

## Verified
* params init (line 139) precedes first use (line 164).
* file compiles (py_compile).

1. We need to creates total of 4 species named A, B, C, D. pymecsim class `Specie` is used to create species with the following information : Diffusion coefficient, initial concentration and whether the specie is expected to be a solution specie.

2. Two reaction we would like to simulate can be defined using the `Reaction` classes. Since we have one `ChargeTransfer` and `ChemicalReaction` we use the corresponding pymecsim classes to build these reactions.

3. Building reactions requires specific information pertaining to the type of reaction. For example charge transfer reaction would need heterogenous reaction constant ks and formal potential E0, Alpha parameter used for a Butler-Volmer model; chemical or catalytic reactions would need forward and backward coeffients. All the reactions by requires dictonaries of reactants and products (see example below for a use case)

4. Once we have set of species and reactions, we can pass them onto a `Mechanism` class to build a mechanism we would like to simulate.

5. We then need to set up our voltage loading for the `Voltammetry` class. pymecsim provides three types of voltage loadings: `DCVoltammetry` that assumes that voltage start and end are same with a sweep reversal at a particualr voltage for a given voltage scan rate. `ACVoltammetry` class that needs to be passed along with a DC voltammetry requires number of AC sources, their amplitude and frequency. Alternatively, user can also pass on a Voltage input file as .txt file (user is referred to the original MECSim documentation for more details)

6. User can also specify an `Electrode` type, which is by default set to be a `PlanarElectrode` with a surface area of 1cm2. Other options include `SphericalElectrode`, `CylindricalElectrode`, `RDEElectrode`.

7. We can then finally create an `Experiment` by passing the following: `Mechanism` along with optional parameters using the keys `electrode` for `Electrode`, `voltammetry` for `Voltammetry` and `capacitance` for `Capacitance`
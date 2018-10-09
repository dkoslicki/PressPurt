# PressPurtCoreAlg
The core algorithm for the Koslicki &amp; Novak JMB paper (2018).

# Installation
I suggest doing this in a [virtual environment](https://docs.python.org/3/library/venv.html).
```bash
pip install -r requirements.txt
```

# Usage
The basic steps are:
1. Preprocess the matrix (`scripts/PreprocessMatrix.py`)
2. Compute statistics on it (`scripts/ComputeEntryWisePerturbationExpectation.py`)
3. Visualize the results (`scripts/GenerateEntryWiseFigures.py`)

# Example
In the `scripts` directory:
```bash
python PreprocessMatrix.py ../ExampleJacobians/IGP.csv ../ExampleJacobians
python ComputeEntryWisePerturbationExpectation.py ../ExampleJacobians
python GenerateEntryWiseFigures.py ../ExampleJacobians -a
```

For quantitative stability:
```bash
python ComputeQuantitativeSensitivity.py ../ExampleJacobians/IGP.csv ../ExampleJacobians/
python GenerateQuantitativeSensitivityFigure.py ../ExampleJacobians/
```

For multi-entry perturbation:
```bash
python ComputeMultiEntryPerturbationExpectation.py ../ExampleJacobians/IGP.csv 
```

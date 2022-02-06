# na2aa
Convert DNA sequences to longest amino acide sequence for all six reading frames

## Installation
I suppose you do not want to keep the package installed, thus I recommend
to install it in a fresh virtual environment.


Clone the repository then change into the directory and execute:
```
pip install .
```

Uninstall:
```
pip uninstall na2aa
```

## Usage

```
na2aa -s <sequences.fasta> -i <intervals.gff>
```

## Example

In the example folder, run:
```
na2aa -s sequences.fasta -i intervals.gff
```
The output is printed to the command line but you can
find the example output also in the file `genes.fasta`.

## Issues

  * make sure the intervals index is 1-based (not 0-based) and follows the same format as in the example

## Future Development
  * add tests with `pytest` 
  * this is just a small showcase project


# Needleman Wunsch Sequence alignment #

Implementation of Needleman-Wunsch sequence alignment algorithm


## Running from terminal ##


### Without matched regions ###

` python3 needleman_wunsch.py -q query_sequence_file -r reference_sequence_file `

The query sequence is human protein sequence, and the reference sequence is the fly protein sequence. Thus, for alignment of HOX protein the following command should be given:

` python3 needleman_wunsch.py -q Human_HOX.fa -r Fly_HOX.fa `


Alignment of PAX protein sequence in human and fruit fly

![picture alt](https://github.com/Mahendra-Maiti/needleman_wunsch/blob/master/align_PAX.png "Alignment of PAX protein")


### Anchored version ###

Anchored alignment of protein sequences in human and fruit fly can be run by adding a file specifying matched regions at the end.

` python3 needleman_wunsch.py -q query_sequence_file -r reference_sequence_file  -m matched_regions`

![picture alt](https://github.com/Mahendra-Maiti/needleman_wunsch/blob/master/align_PAX.png "Alignment of PAX protein")


# OpenMP_word
## Task Description:
- Each node reads a portion of the text and the dictionary file, with the goal of finding misspelled words in the text. The misspelled words differ by only one character from the correct words in the dictionary.
- Nodes use MergeSort to sort the text for efficient comparison.
- After each node finishes comparing its assigned dictionary with the text, it passes the dictionary to the next node in the sequence.
- This process continues until each node has compared its text portion with all parts of the dictionary, ensuring a comprehensive spell-check across all nodes.
output Example:
appel: apple


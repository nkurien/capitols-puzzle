# The Capitols puzzle

The 11-by-11 region below can be subdivided into various states. A "state" is a connected group of 1 or more squares that is symmetric in some way, be it via rotation or reflection. If every nontrivial rotation or reflection of the plane that takes a state to itself leaves exactly 1 of its squares fixed, then that square is its capitol. (So, for example, a 1-by-3 region has no capitol, nor does the "X" pentomino.)

The cost of traveling from one square to an orthogonally adjacent square in this nation is equal to the size (i.e. number of squares) of the state you are travelling in, or, if travelling between two states, the lesser of the two states' sizes. A number in a cell (below) denotes the lowest cost of travelling from the numbered square to any state capitol.

After determining the shapes of the states, label each cell with the size of the state it is in. Compute the sum of these values in each row, and submit as your answer the sum of the squares of these row sums (as in the example).

## The puzzle (`python solve.py main`)

|   |   |   |   |   |   |   |   |   |   |   |
|---|---|---|---|---|---|---|---|---|---|---|
|   | 8 | 4 |   |   | 11 |   |   |   | 5 |   |
|   |   |   |   |   |   |   | 11 |   |   |   |
|   |   | 7 |   |   | 1 |   |   |   |   | 14 |
| 28 |   |   | 51 |   |   | 1 |   | 6 |   |   |
|   | 22 |   |   |   |   |   |   | 4 |   | 1 |
|   |   |   | 15 |   | 10 |   | 0 |   |   |   |
| 11 |   | 11 |   |   |   |   |   |   | 9 |   |
|   |   | 17 |   | 14 |   |   | 10 |   |   | 13 |
| 30 |   |   |   |   | 6 |   |   | 45 |   |   |
|   |   |   | 10 |   |   |   |   |   |   |   |
|   | 26 |   |   |   | 0 |   |   | 77 | 61 |   |

## Worked example (`python solve.py ex`)

Example grid:

|   |   |   |   |   |
|---|---|---|---|---|
| 12 |   |   |   | 17 |
| 8 |   | 6 | 3 |   |
|   |   |   |   |   |
|   | 7 | 6 |   | 10 |
| 25 |   |   |   | 0 |

Solved, with each cell labelled by the size of its state (letters mark the states, `*` marks a capitol):

```
A B B C C        7  4  4 11 11
A B B D C        7  4  4  3 11
A A* D D* C      7  7  3  3 11
C A A A C        11 7  7  7 11
C C C C C*       11 11 11 11 11
```

Row sums: 37, 29, 31, 43, 55. Answer: 37² + 29² + 31² + 43² + 55² = **8,045**.

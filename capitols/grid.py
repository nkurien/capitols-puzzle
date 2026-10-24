class Grid:
    """An N x N board with cells indexed row-major."""

    def __init__(self, N):
        self.N = N
        self.cells = [(r, c) for r in range(N) for c in range(N)]
        self.idx = {x: i for i, x in enumerate(self.cells)}
        self.n = len(self.cells)
        self.adj = [list(self._neighbours(i)) for i in range(self.n)]

    def _neighbours(self, i):
        r, c = self.cells[i]
        for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            if (r + dr, c + dc) in self.idx:
                yield self.idx[(r + dr, c + dc)]

    def edges(self):
        return [(i, j) for i in range(self.n) for j in self.adj[i] if i < j]


def monotone_paths(c1, c2):
    """Every shortest orthogonal path from c1 to c2, as lists of (from, to) steps."""
    r2, c2_col = c2
    paths = []

    def dfs(curr, path):
        if curr == c2:
            paths.append(path)
            return
        cr, cc = curr
        if cr != r2:
            dr = 1 if r2 > cr else -1
            nxt = (cr + dr, cc)
            dfs(nxt, path + [(curr, nxt)])
        if cc != c2_col:
            dc = 1 if c2_col > cc else -1
            nxt = (cr, cc + dc)
            dfs(nxt, path + [(curr, nxt)])

    dfs(c1, [])
    return paths

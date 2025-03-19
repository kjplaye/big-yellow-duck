import graphviz
import numpy as np
import tempfile
from cat import Cat
from pca import PCA

class Graph():
    """A directed graph with hooks.
         
    edges is a list of pairs of hashables
    
    USAGE EXAMPLE:
    >>> G = Graph([[10,11],[10,12],[11,12],[12,11]])
    >>> print(G)
    >>> print("(edge, vertex) adjacency:")
    >>> print(G.matrix())
    >>> G.dot(color = {10:(1.0,0.7,0.0)},
            fill_color = {12:(0.5,0.2,0.9)},
            edge_color = {0:(0.2,0.8,0.2)},
            edge_label = {1:"bark"},
            edge_size = {1:4.0})
    """
    def __init__(self, edges):
        self._edge_cat = Cat([x for e in edges for x in e])
        self.vertices = self._edge_cat.cat
        self.edges = [self._edge_cat[i:i+2] for i in range(0,len(self._edge_cat),2)]
    def matrix(self, return_type = 'ev', directed = True):
        """Return adjacency matrix.

        G.matrix()              is (edge,   vertex) adjacency (return_type = 'ev')
        G.matrix().T            is (vertex,   edge) adjacency (return_type = 've')
        G.matrix().T @ G.matrix() is (vertex, vertex) adjacency (return_type = 'vv')

        `directed = False` changes edge rows from [1,-1] to [1,1].  The (vertex, vertex) flips
        sign off the diagonal.
        """
        M = np.zeros([len(self.edges),len(self.vertices)])
        for i,e in enumerate(self.edges):
            M[i,e[0]] = 1
            M[i,e[1]] = -1 if directed else 1
        if return_type == 'ev':
            return M
        elif return_type == 've':
            return M.T
        elif return_type == 'vv':
            return M.T @ M
        raise NotImplementedError("return_type not expected")
    def spectral(self, matrix_return_type = 'ev', directed = True):
        """Use PCA on Incidence Matrix."""
        return PCA(self.matrix(return_type = matrix_return_type, directed = directed))
    def dot(self, color = {}, fill_color = {}, edge_color = {}, edge_label = {},
            edge_size = {}, comment = 'Dot from Python', ranksep = None):
        def _color_code(color):
            return "#%02x%02x%02x" % tuple((int(255 * color[i]) for i in range(3)))                      
        GV = graphviz.Digraph(comment=comment)
        if ranksep is not None:
            GV.graph_attr['ranksep'] = str(ranksep)
        for v in self.vertices:
            c = "#000000" if v not in color else _color_code(color[v])
            fc = "#ffffff" if v not in fill_color else _color_code(fill_color[v])
            GV.node(str(v), color=c, fillcolor = fc, style= 'filled')
        for k, e in enumerate(self.edges):
            ec = "#000000" if k not in edge_color else _color_code(edge_color[k])
            es = "1.0" if k not in edge_size else "%f" % edge_size[k]
            el = "" if k not in edge_label else edge_label[k]
            GV.edge(str(self.vertices[e[0]]),str(self.vertices[e[1]]),
                    color = ec, penwidth = es, label = el)
        f = tempfile.NamedTemporaryFile(suffix = f'_{comment}.pdf')
        GV.render(f.name[:-4], view=True)
        print("dot temp file:", f.name)
        input("Press key to continue...")
    def count_cycles(self, directed = True):
        if not(directed):
            return len(self.edges) - len(self.vertices) + len(self.connected_components())
        num_nodes = len(self.vertices)
        def dfs(node, visited, path):
            # Set some things
            visited[node] = True
            path[node] = True
            count = 0
            # Recurse into all next node to extend path
            for neighbor in graph[node]:
                if not visited[neighbor]:
                    count += dfs(neighbor, visited, path)
                elif path[neighbor]:  # Cycle detected
                    count += 1
                # Back up path
                path[node] = False
            return count
        # Set up stuff
        num_cycles = 0
        visited = [False] * num_nodes
        path = [False] * num_nodes    
        graph = [[] for _ in range(num_nodes)]
        for edge in self.edges:
            graph[edge[0]].append(edge[1])
        # Start recursion
        for node in range(num_nodes):
            if not visited[node]:
                num_cycles += dfs(node, visited, path)    
        return num_cycles
    def connected_components(self):
        visited = set()
        components = []
        G = {v:[] for v in set(self._edge_cat)}
        for edge in self.edges:
            G[edge[0]].append(edge[1])
        for node in G:
            if node not in visited:
                component = set()
                queue = [node]    
                while queue:
                    current_node = queue.pop(0)
                    if current_node not in visited:
                        component.add(current_node)
                        visited.add(current_node)
                        queue.extend(G.get(current_node, []))    
                    components.append(component)    
        return set(map(tuple, components))       
    def __repr__(self):
        return f'Directed graph: {len(self.vertices)} vertices, {len(self.edges)} edges'

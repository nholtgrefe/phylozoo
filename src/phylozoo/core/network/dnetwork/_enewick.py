"""
eNewick (Extended Newick) format parser and writer.

This module provides complete eNewick functionality:
- Parsing eNewick strings into a generic intermediate representation
- Converting DirectedPhyNetwork objects to eNewick format strings
- Converting eNewick strings to DirectedPhyNetwork objects

eNewick format extends Newick format to support phylogenetic networks with
hybrid nodes. Hybrid nodes are marked with #H1, #H2, etc. in the string.

eNewick Format Specification
----------------------------

The eNewick format extends standard Newick format with the following features:

1. **Tree Structure**: Uses parentheses to represent nested tree structure
   - `(A,B)` represents a node with two children A and B
   - `((A,B),C)` represents a nested structure

2. **Node Labels**: 
   - Unquoted: `A`, `Species1`, `node_123` (alphanumeric, underscore, dash, dot)
   - Quoted: `'Species Name'`, `'Node with spaces'` (single quotes, escape quotes with '')
   - Internal nodes can be unlabeled (auto-generated IDs)

3. **Branch Lengths**: Specified after nodes with `:value`
   - `A:0.5` - edge to A has length 0.5
   - `(A:0.5,B:0.3):0.1` - edges to A and B have lengths 0.5 and 0.3, 
     edge to parent has length 0.1
   - Supports scientific notation: `A:1.5e-3`

4. **Hybrid Nodes**: Marked with `#H1`, `#H2`, etc. (Extended Newick format)
   - `(A,B)#H1` - defines hybrid node #1 with children A and B
   - `#H1` - references the previously defined hybrid node #1
   - A hybrid node can appear multiple times (once with definition, then one or more references)
   - Each reference creates an additional parent edge to the hybrid node (reticulation)
   - Example: `((A,B)#H1,#H1,#H1);` creates a hybrid with 3 parents

5. **Comments**: Enclosed in square brackets `[comment]`
   - `A[comment]` - comment attached to node A (after label)
   - `[comment]A` - comment before node label
   - Comments are stored as 'comment' attribute in node dictionaries
   - Note: Comments can appear before or after labels, but not after branch lengths
     or hybrid markers in the current implementation

6. **Termination**: Must end with semicolon `;`

Examples
--------

Basic tree:
    >>> result = parse_enewick("((A,B),C);")
    >>> len(result.edges)
    4
    >>> result.root  # Internal node IDs are integers
    0

Tree with branch lengths:
    >>> result = parse_enewick("((A:0.5,B:0.3):0.1,C:0.2);")
    >>> result.edges[0]['branch_length']
    0.5
    >>> result.edges[2]['branch_length']  # Edge to internal node
    0.1

Network with hybrid node (Extended Newick format):
    >>> result = parse_enewick("((A,B)#H1,#H1);")
    >>> result.hybrid_nodes
    {1: 1}
    >>> # Hybrid node appears twice: once with definition, once as reference
    >>> # This creates 2 incoming edges to the hybrid node

Quoted node labels:
    >>> result = parse_enewick("('Species A','Species B');")
    >>> [n['label'] for n in result.nodes if 'label' in n]
    ['Species A', 'Species B']

Comments (before or after labels):
    >>> result = parse_enewick("(A[comment1],B[comment2]);")
    >>> result.nodes[0].get('comment')
    'comment1'
    >>> result = parse_enewick("([comment1]A,[comment2]B);")
    >>> result.nodes[0].get('comment')
    'comment1'

Complex example with multiple features:
    >>> enewick = "((('Leaf A':0.5,'Leaf B':0.3)#H1:0.1),'Leaf C':0.2);"
    >>> result = parse_enewick(enewick)
    >>> # Contains:
    >>> # - Quoted labels: 'Leaf A', 'Leaf B', 'Leaf C'
    >>> # - Branch lengths on edges: 0.5, 0.3, 0.1, 0.2
    >>> # - Hybrid node: #H1 (stored in result.hybrid_nodes)
    >>> len([e for e in result.edges if 'branch_length' in e])
    4
    >>> result.hybrid_nodes
    {'internal_1': 1}

Scientific notation for branch lengths:
    >>> result = parse_enewick("(A:1.5e-3,B:2.0e2);")
    >>> result.edges[0]['branch_length']
    0.0015
    >>> result.edges[1]['branch_length']
    200.0

Unlabeled internal nodes (auto-generated IDs):
    >>> result = parse_enewick("((A,B),C);")
    >>> [n['id'] for n in result.nodes if isinstance(n['id'], int)]
    [0, 1]

Non-binary nodes (polytomies):
    >>> result = parse_enewick("(A,B,C,D);")
    >>> # Root has 4 children
    >>> len([e for e in result.edges if e['u'] == result.root])
    4

Hybrid node with multiple parents:
    >>> result = parse_enewick("((A,B)#H1,#H1,#H1);")
    >>> # Hybrid appears 3 times (1 definition + 2 references) = 3 parents
    >>> hybrid_id = list(result.hybrid_nodes.keys())[0]
    >>> len([e for e in result.edges if e['v'] == hybrid_id])
    3

Notes
-----
- Internal nodes without labels are assigned auto-generated integer IDs (0, 1, 2, ...)
- Leaves use their label as both `id` and `label` (strings)
- Branch lengths are stored on edges, not nodes
- Hybrid numbers (#H1, #H2, etc.) are tracked in the `hybrid_nodes` dictionary
- Comments are stored as 'comment' attribute in node dictionaries
- The parser is case-sensitive for node labels
- Non-binary nodes (polytomies) are fully supported: `(A,B,C,D)` creates a node with 4 children
- Hybrid nodes can have multiple parents by appearing multiple times: `((A)#H1,#H1,#H1)` creates 3 parent edges
- All labels must be unique; duplicate labels raise ``ENewickParseError``
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any


@dataclass
class ParsedENewick:
    """
    Intermediate representation of a parsed eNewick string.
    
    This structure is format-agnostic and can be converted to any network type
    (DirectedPhyNetwork, SemiDirectedPhyNetwork, etc.).
    
    Attributes
    ----------
    edges : list[dict[str, Any]]
        List of edges, each as a dict with keys:
        - 'u': source node ID
        - 'v': target node ID
        - 'key': edge key (for parallel edges, default 0)
        - Additional edge attributes (branch_length, bootstrap, gamma, etc.)
    nodes : list[dict[str, Any]]
        List of nodes, each as a dict with keys:
        - 'id': node ID
        - 'label': node label (if present)
        - Additional node attributes
    root : Any | None
        Root node ID, or None if not determined
    hybrid_nodes : dict[Any, int]
        Mapping from hybrid node IDs to their hybrid numbers (#H1, #H2, etc.)
    metadata : dict[str, Any]
        Additional metadata from the eNewick string
    """
    edges: list[dict[str, Any]] = field(default_factory=list)
    nodes: list[dict[str, Any]] = field(default_factory=list)
    root: Any | None = None
    hybrid_nodes: dict[Any, int] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)


class ENewickParseError(Exception):
    """Exception raised when eNewick parsing fails."""
    pass


def parse_enewick(enewick_string: str) -> ParsedENewick:
    """
    Parse an eNewick string into a generic intermediate representation.
    
    See the module docstring for detailed format specification and examples.
    
    Parameters
    ----------
    enewick_string : str
        The eNewick string to parse. Must end with ';'.
    
    Returns
    -------
    ParsedENewick
        Parsed representation with edges, nodes, root, and hybrid node information.
        - `edges`: List of edge dictionaries with 'u', 'v', 'key', and optional
          attributes like 'branch_length'
        - `nodes`: List of node dictionaries with 'id', optional 'label', and
          optional attributes like 'comment'
        - `root`: Root node ID
        - `hybrid_nodes`: Dictionary mapping hybrid node IDs to hybrid numbers
        - `metadata`: Additional metadata dictionary
    
    Raises
    ------
    ENewickParseError
        If the eNewick string is malformed or cannot be parsed.
    
    Examples
    --------
    >>> # Simple tree
    >>> result = parse_enewick("((A,B),C);")
    >>> len(result.edges)
    4
    >>> result.root
    0
    
    >>> # Tree with branch lengths
    >>> result = parse_enewick("((A:0.5,B:0.3):0.1,C:0.2);")
    >>> result.edges[0]['branch_length']
    0.5
    
    >>> # Network with hybrid node
    >>> result = parse_enewick("((A,B)#H1,C);")
    >>> result.hybrid_nodes
    {1: 1}
    
    >>> # Quoted labels and comments
    >>> result = parse_enewick("('Species A'[note1],'Species B'[note2]);")
    >>> [n['label'] for n in result.nodes if 'label' in n]
    ['Species A', 'Species B']
    >>> result.nodes[0].get('comment')
    'note1'
    """
    # Remove whitespace and check for semicolon
    enewick_string = enewick_string.strip()
    if not enewick_string.endswith(';'):
        raise ENewickParseError("eNewick string must end with ';'")
    
    enewick_string = enewick_string[:-1].strip()  # Remove semicolon
    
    # Initialize parser state
    parser = _ENewickParser(enewick_string)
    root_id, _ = parser.parse()
    
    # Update result with root
    result = ParsedENewick(
        edges=parser.edges,
        nodes=parser.nodes,
        root=root_id,
        hybrid_nodes=parser.hybrid_nodes,
        metadata={}
    )
    
    return result


class _ENewickParser:
    """
    Internal parser class for eNewick strings.
    
    Uses recursive descent parsing to handle the nested structure.
    """
    
    def __init__(self, enewick_string: str):
        self.enewick_string = enewick_string
        self.pos = 0
        self.length = len(enewick_string)
        
        # Track nodes and edges as we parse
        self.edges: list[dict[str, Any]] = []
        self.nodes: list[dict[str, Any]] = []
        self.hybrid_nodes: dict[Any, int] = {}
        self.node_counter = 0  # For generating internal node IDs (integers)
        self.hybrid_counter = 0  # For tracking hybrid numbers
        self._hybrid_id_map: dict[int, Any] = {}  # hybrid number -> node id
        self._labels_seen: set[str] = set()  # Track labels to prevent duplicates
        
        # Track parent-child relationships
        self.parent_stack: list[Any] = []
        
    def parse(self) -> tuple[Any, float | None]:
        """
        Parse the eNewick string.
        
        Returns
        -------
        tuple[Any, float | None]
            (root_id, root_branch_length) - root has no parent so branch_length is None.
        """
        root_id, root_branch_length = self._parse_subtree()
        
        if self.pos < self.length:
            raise ENewickParseError(
                f"Unexpected characters after tree at position {self.pos}: "
                f"'{self.enewick_string[self.pos:]}'"
            )
        
        return root_id, root_branch_length
    
    def _parse_subtree(self) -> tuple[Any, float | None]:
        """
        Parse a subtree (node with optional children).
        
        Returns
        -------
        tuple[Any, float | None]
            (node_id, branch_length) where branch_length is for edge TO this node from parent.
        """
        # Check if this is a leaf or internal node
        if self._peek() == '(':
            # Internal node: parse children
            return self._parse_internal_node()
        else:
            # Leaf node
            return self._parse_leaf()
    
    def _parse_internal_node(self) -> tuple[Any, float | None]:
        """
        Parse an internal node with children.
        
        Returns
        -------
        tuple[Any, float | None]
            (node_id, branch_length) where branch_length is for edge TO this node from parent.
        """
        self._skip_whitespace()
        # Consume opening parenthesis
        self._expect('(')
        self._skip_whitespace()
        
        # Generate internal node ID (integer, always unique)
        internal_id = self.node_counter
        self.node_counter += 1
        
        # Parse children (each child subtree returns (child_id, branch_length))
        children: list[tuple[Any, float | None]] = []
        first_child = True
        
        while self._peek() != ')':
            self._skip_whitespace()
            if not first_child:
                self._expect(',')
                self._skip_whitespace()
            
            # Parse child subtree (returns child_id and its branch_length)
            child_id, child_branch_length = self._parse_subtree()
            
            children.append((child_id, child_branch_length))
            first_child = False
        
        # Consume closing parenthesis
        self._expect(')')
        self._skip_whitespace()
        
        # Parse node label and attributes
        node_label, node_attrs = self._parse_node_label_and_attrs()
        
        # Parse branch length (if present) - this is for edge TO this node from parent
        parent_branch_length = self._parse_branch_length()
        
        # Parse hybrid marker (if present)
        hybrid_number = self._parse_hybrid_marker()
        
        # Handle hybrid nodes: reuse existing node ID if hybrid number seen before
        is_hybrid_duplicate = False
        if hybrid_number is not None:
            if hybrid_number in self._hybrid_id_map:
                internal_id = self._hybrid_id_map[hybrid_number]
                is_hybrid_duplicate = True
            else:
                self._hybrid_id_map[hybrid_number] = internal_id

        # Enforce unique labels (for provided labels only)
        if node_label:
            if node_label in self._labels_seen:
                raise ENewickParseError(f"Duplicate label detected: '{node_label}'")
            self._labels_seen.add(node_label)

        # Record node (internal_id stays integer; label stored separately)
        # For duplicate hybrid occurrences, avoid duplicating the node entry
        if not is_hybrid_duplicate:
            node_data: dict[str, Any] = {'id': internal_id}
            if node_label:
                node_data['label'] = node_label
            node_data.update(node_attrs)
            self.nodes.append(node_data)
        else:
            # If duplicate hybrid occurrence has a label/attrs, ensure consistency
            if node_label:
                # If label differs from an existing one, raise error
                existing = next((n for n in self.nodes if n.get('id') == internal_id), None)
                if existing is not None and existing.get('label') and existing.get('label') != node_label:
                    raise ENewickParseError(
                        f"Hybrid node #{hybrid_number} label mismatch: "
                        f"'{existing.get('label')}' vs '{node_label}'"
                    )

        # Record hybrid node
        if hybrid_number is not None:
            self.hybrid_nodes[internal_id] = hybrid_number
        
        # Create edges from internal node to children
        # Branch length on edge is the branch_length returned by the child subtree
        for child_id, child_branch_length in children:
            edge_data: dict[str, Any] = {
                'u': internal_id,
                'v': child_id,
                'key': 0
            }
            if child_branch_length is not None:
                edge_data['branch_length'] = child_branch_length
            self.edges.append(edge_data)
        
        return internal_id, parent_branch_length
    
    def _parse_leaf(self) -> tuple[Any, float | None]:
        """
        Parse a leaf node.
        
        Returns
        -------
        tuple[Any, float | None]
            (node_id, branch_length) where branch_length is for edge TO this leaf from parent.
        """
        self._skip_whitespace()
        
        # Check if this is a bare hybrid node reference (starts with #H)
        if self._peek() == '#':
            # This is a reference to an existing hybrid node
            hybrid_number = self._parse_hybrid_marker()
            if hybrid_number is None:
                raise ENewickParseError("Expected hybrid marker after '#'")
            
            # Look up the existing hybrid node
            if hybrid_number not in self._hybrid_id_map:
                raise ENewickParseError(
                    f"Hybrid node reference #H{hybrid_number} found before definition"
                )
            
            node_id = self._hybrid_id_map[hybrid_number]
            
            # Parse branch length (if present)
            branch_length = self._parse_branch_length()
            
            # Don't create a new node entry - just return the existing node ID
            return node_id, branch_length
        
        # Parse node label
        node_label, node_attrs = self._parse_node_label_and_attrs()
        
        if not node_label:
            # Generate leaf ID/label if no label provided (string)
            node_label = f"leaf_{self.node_counter}"
            self.node_counter += 1

        # Enforce unique labels for leaves
        if node_label in self._labels_seen:
            raise ENewickParseError(f"Duplicate label detected: '{node_label}'")
        self._labels_seen.add(node_label)
        
        # Parse branch length (if present) - this is for edge TO this leaf from parent
        branch_length = self._parse_branch_length()
        
        # Parse hybrid marker (if present) - leaves can be hybrid nodes
        hybrid_number = self._parse_hybrid_marker()
        
        # Record node (leaves use label as both id and label)
        node_data: dict[str, Any] = {'id': node_label, 'label': node_label}
        node_data.update(node_attrs)
        self.nodes.append(node_data)
        
        # Record hybrid node (first occurrence - define it)
        if hybrid_number is not None:
            if hybrid_number in self._hybrid_id_map:
                raise ENewickParseError(
                    f"Hybrid node #H{hybrid_number} defined multiple times with children"
                )
            self.hybrid_nodes[node_label] = hybrid_number
            self._hybrid_id_map[hybrid_number] = node_label
        
        return node_label, branch_length
    
    def _parse_node_label_and_attrs(self) -> tuple[str | None, dict[str, Any]]:
        """
        Parse node label and attributes.
        
        Returns
        -------
        tuple[str | None, dict[str, Any]]
            (label, attributes_dict)
        """
        self._skip_whitespace()
        attrs: dict[str, Any] = {}
        label: str | None = None
        
        # Check for comment (square brackets) before label
        if self._peek() == '[':
            comment = self._parse_comment()
            attrs['comment'] = comment
            self._skip_whitespace()
        
        # Parse label (identifier or quoted string)
        if self._peek() not in (':', ')', ',', ';', '#', '['):
            label = self._parse_identifier()
            self._skip_whitespace()
        
        # Check for attributes in square brackets after label
        if self._peek() == '[':
            comment = self._parse_comment()
            attrs['comment'] = comment
            self._skip_whitespace()
        
        return label, attrs
    
    def _parse_identifier(self) -> str:
        """
        Parse an identifier (node label).
        
        Can be unquoted (alphanumeric + underscore) or quoted string.
        """
        if self._peek() == "'":
            return self._parse_quoted_string()
        
        # Unquoted identifier: alphanumeric, underscore, dash, dot
        # Stop at special characters: : ) , ; # [
        start = self.pos
        while (self.pos < self.length and 
               self.enewick_string[self.pos] not in (':', ')', ',', ';', '#', '[', '(', ' ')):
            self.pos += 1
        
        if self.pos == start:
            return ""
        
        return self.enewick_string[start:self.pos].strip()
    
    def _parse_quoted_string(self) -> str:
        """Parse a quoted string (single quotes)."""
        self._expect("'")
        start = self.pos
        
        # Find closing quote (handle escaped quotes)
        while self.pos < self.length:
            if self.enewick_string[self.pos] == "'":
                if self.pos + 1 < self.length and self.enewick_string[self.pos + 1] == "'":
                    # Escaped quote
                    self.pos += 2
                else:
                    # Closing quote
                    result = self.enewick_string[start:self.pos]
                    self.pos += 1
                    return result.replace("''", "'")
            self.pos += 1
        
        raise ENewickParseError("Unclosed quoted string")
    
    def _parse_branch_length(self) -> float | None:
        """Parse branch length (e.g., :0.5)."""
        self._skip_whitespace()
        if self._peek() != ':':
            return None
        
        self.pos += 1  # Consume ':'
        self._skip_whitespace()
        
        # Parse number (can be negative, scientific notation, etc.)
        start = self.pos
        
        # Optional sign
        if self.pos < self.length and self.enewick_string[self.pos] in ('-', '+'):
            self.pos += 1
        
        # Parse digits and decimal point
        has_digit = False
        while (self.pos < self.length and 
               (self.enewick_string[self.pos].isdigit() or self.enewick_string[self.pos] == '.')):
            if self.enewick_string[self.pos].isdigit():
                has_digit = True
            self.pos += 1
        
        # Optional scientific notation
        if (self.pos < self.length and 
            self.enewick_string[self.pos] in ('e', 'E') and has_digit):
            self.pos += 1
            if self.pos < self.length and self.enewick_string[self.pos] in ('-', '+'):
                self.pos += 1
            while self.pos < self.length and self.enewick_string[self.pos].isdigit():
                self.pos += 1
        
        if self.pos == start or not has_digit:
            raise ENewickParseError(f"Expected number after ':' at position {start}")
        
        try:
            return float(self.enewick_string[start:self.pos])
        except ValueError:
            raise ENewickParseError(
                f"Invalid branch length at position {start}: "
                f"'{self.enewick_string[start:self.pos]}'"
            )
    
    def _parse_hybrid_marker(self) -> int | None:
        """Parse hybrid marker (e.g., #H1)."""
        self._skip_whitespace()
        if self._peek() != '#':
            return None
        
        self.pos += 1  # Consume '#'
        
        # Expect 'H' followed by number
        if self.pos >= self.length or self.enewick_string[self.pos] != 'H':
            raise ENewickParseError(
                f"Expected 'H' after '#' at position {self.pos}"
            )
        
        self.pos += 1  # Consume 'H'
        
        # Parse hybrid number
        start = self.pos
        while self.pos < self.length and self.enewick_string[self.pos].isdigit():
            self.pos += 1
        
        if self.pos == start:
            raise ENewickParseError(
                f"Expected number after '#H' at position {self.pos}"
            )
        
        hybrid_num = int(self.enewick_string[start:self.pos])
        
        # Track hybrid numbers
        if hybrid_num not in self.hybrid_nodes.values():
            self.hybrid_counter = max(self.hybrid_counter, hybrid_num)
        
        return hybrid_num
    
    def _parse_comment(self) -> str:
        """Parse comment in square brackets [comment]."""
        self._expect('[')
        start = self.pos
        
        # Find closing bracket
        depth = 1
        while self.pos < self.length and depth > 0:
            if self.enewick_string[self.pos] == '[':
                depth += 1
            elif self.enewick_string[self.pos] == ']':
                depth -= 1
            self.pos += 1
        
        if depth > 0:
            raise ENewickParseError("Unclosed comment bracket")
        
        # Remove closing bracket from result
        return self.enewick_string[start:self.pos - 1]
    
    def _peek(self) -> str | None:
        """Peek at next character without consuming it."""
        if self.pos >= self.length:
            return None
        return self.enewick_string[self.pos]
    
    def _expect(self, char: str) -> None:
        """Expect and consume a specific character."""
        if self.pos >= self.length:
            raise ENewickParseError(
                f"Unexpected end of string, expected '{char}'"
            )
        if self.enewick_string[self.pos] != char:
            raise ENewickParseError(
                f"Expected '{char}' at position {self.pos}, "
                f"got '{self.enewick_string[self.pos]}'"
            )
        self.pos += 1
    
    def _skip_whitespace(self) -> None:
        """Skip whitespace characters."""
        while self.pos < self.length and self.enewick_string[self.pos].isspace():
            self.pos += 1



# ========== eNewick Writer and Reader for DirectedPhyNetwork ==========

from typing import TYPE_CHECKING, Any, TypeVar

if TYPE_CHECKING:
    from .base import DirectedPhyNetwork

T = TypeVar('T')


def to_enewick(network: 'DirectedPhyNetwork', **kwargs: Any) -> str:
    """
    Convert a DirectedPhyNetwork to Extended Newick (eNewick) format.
    
    This function serializes a directed phylogenetic network to the Extended Newick
    format, which supports hybrid nodes (reticulations) using the #H marker notation.
    
    Features:
    - Branch lengths are encoded as :length (e.g., :0.5)
    - Gamma and bootstrap values are encoded as comments: [&gamma=0.6,bootstrap=0.95]
    - Internal node labels are included when present
    - Hybrid nodes use Extended Newick #Hn markers
    - Output is deterministic (children sorted by node ID then label)
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to serialize.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    str
        The eNewick string representation of the network.
    
    Raises
    ------
    ValueError
        If the network is empty (has no nodes).
    
    Examples
    --------
    >>> # Simple tree
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> to_enewick(net)
    '(A,B);'
    
    >>> # Tree with branch lengths
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         {'u': 3, 'v': 1, 'branch_length': 0.5},
    ...         {'u': 3, 'v': 2, 'branch_length': 0.3}
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> to_enewick(net)
    '(A:0.5,B:0.3);'
    
    Notes
    -----
    - For parallel edges between the same nodes, only the first edge is used
    - Labels containing special characters (spaces, parentheses, colons, etc.) 
      are automatically quoted with single quotes
    - The output is deterministic: multiple calls on the same network produce 
      the same string
    """
    if network.number_of_nodes() == 0:
        raise ValueError("Cannot convert empty network to eNewick")
    
    # Single node network
    if network.number_of_nodes() == 1:
        node = next(iter(network._graph.nodes))
        label = network.get_label(node)
        if label:
            return f"{_quote_label_if_needed(label)};"
        else:
            return f"{node};"
    
    # Initialize hybrid tracking
    hybrid_nodes_set = set(network.hybrid_nodes)
    hybrid_to_id: dict[T, int] = {}
    hybrid_counter = 1
    
    # Assign IDs to hybrid nodes deterministically
    for hybrid in sorted(hybrid_nodes_set, key=lambda n: (str(n), network.get_label(n) or '')):
        hybrid_to_id[hybrid] = hybrid_counter
        hybrid_counter += 1
    
    # Track which hybrids have been defined
    defined_hybrids: set[T] = set()
    
    def build_subtree(node: T, parent_edge_data: dict[str, Any] | None = None) -> str:
        """
        Recursively build eNewick string for subtree rooted at node.
        
        Parameters
        ----------
        node : T
            Current node to process.
        parent_edge_data : dict[str, Any] | None
            Edge data from parent to this node (for branch length, gamma, bootstrap).
        
        Returns
        -------
        str
            eNewick substring for this subtree.
        """
        # Check if this is a hybrid node
        is_hybrid = node in hybrid_nodes_set
        
        # If hybrid and already defined, return reference only
        if is_hybrid and node in defined_hybrids:
            # Just return the reference #Hn with parent edge attributes
            hybrid_id = hybrid_to_id[node]
            result = f"#H{hybrid_id}"
            result += _format_edge_attributes(parent_edge_data)
            return result
        
        # Mark hybrid as defined if this is first occurrence
        if is_hybrid:
            defined_hybrids.add(node)
        
        # Get children (sorted deterministically)
        children_list = list(network.children(node))
        children_sorted = sorted(children_list, key=lambda n: (str(n), network.get_label(n) or ''))
        
        # If this is a leaf, just return the label
        if len(children_sorted) == 0:
            label = network.get_label(node)
            if label:
                result = _quote_label_if_needed(label)
            else:
                result = str(node)
            
            # Add hybrid marker if this is a hybrid leaf
            if is_hybrid:
                hybrid_id = hybrid_to_id[node]
                result += f"#H{hybrid_id}"
            
            # Add parent edge attributes
            result += _format_edge_attributes(parent_edge_data)
            return result
        
        # Internal node: process children
        child_strings = []
        for child in children_sorted:
            # Get edge data from node to child
            # Handle parallel edges by using the first edge
            edge_data = _get_first_edge_data(network, node, child)
            child_str = build_subtree(child, edge_data)
            child_strings.append(child_str)
        
        # Build the result string
        result = f"({','.join(child_strings)})"
        
        # Add internal node label if present
        label = network.get_label(node)
        if label:
            result += _quote_label_if_needed(label)
        
        # Add hybrid marker if this is a hybrid
        if is_hybrid:
            hybrid_id = hybrid_to_id[node]
            result += f"#H{hybrid_id}"
        
        # Add parent edge attributes
        result += _format_edge_attributes(parent_edge_data)
        
        return result
    
    # Start from root
    root = network.root_node
    enewick_str = build_subtree(root, None)
    
    return enewick_str + ";"


def _get_first_edge_data(network: 'DirectedPhyNetwork', u: T, v: T) -> dict[str, Any]:
    """
    Get edge data for the first edge from u to v.
    
    For parallel edges, returns data from the first edge encountered.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The network.
    u, v : T
        Edge endpoints.
    
    Returns
    -------
    dict[str, Any]
        Edge data dictionary.
    """
    # Use incident_child_edges to get all edges from u
    for edge_tuple in network.incident_child_edges(u, keys=True, data=True):
        if len(edge_tuple) == 4:
            edge_u, edge_v, key, data = edge_tuple
            if edge_v == v:
                return dict(data) if data else {}
    
    return {}


def _format_edge_attributes(edge_data: dict[str, Any] | None) -> str:
    """
    Format edge attributes (branch_length, gamma, bootstrap) for eNewick.
    
    Parameters
    ----------
    edge_data : dict[str, Any] | None
        Edge data dictionary.
    
    Returns
    -------
    str
        Formatted attribute string (e.g., "[&gamma=0.6,bootstrap=0.95]:0.5").
    """
    if edge_data is None:
        return ""
    
    result = ""
    
    # Build comment section for gamma and bootstrap
    comment_parts = []
    if 'gamma' in edge_data:
        gamma_val = edge_data['gamma']
        comment_parts.append(f"gamma={gamma_val}")
    if 'bootstrap' in edge_data:
        bootstrap_val = edge_data['bootstrap']
        comment_parts.append(f"bootstrap={bootstrap_val}")
    
    if comment_parts:
        result += f"[&{','.join(comment_parts)}]"
    
    # Add branch length
    if 'branch_length' in edge_data:
        branch_length = edge_data['branch_length']
        result += f":{branch_length}"
    
    return result


def _quote_label_if_needed(label: str) -> str:
    """
    Quote a label if it contains special characters.
    
    Special characters that require quoting: spaces, parentheses, colons,
    semicolons, commas, brackets, quotes.
    
    Parameters
    ----------
    label : str
        The label to potentially quote.
    
    Returns
    -------
    str
        The label, quoted if necessary.
    """
    special_chars = {' ', '(', ')', ':', ';', ',', '[', ']', "'", '"', '#'}
    
    if any(char in label for char in special_chars):
        # Escape single quotes in the label
        escaped_label = label.replace("'", "''")
        return f"'{escaped_label}'"
    
    return label


def _parse_comment_attributes(comment: str) -> dict[str, Any]:
    """
    Parse comment string to extract edge attributes (gamma, bootstrap).
    
    Comments in eNewick format can contain edge attributes in the format:
    [&gamma=0.6,bootstrap=0.95]
    
    Parameters
    ----------
    comment : str
        Comment string, possibly starting with '&'.
    
    Returns
    -------
    dict[str, Any]
        Dictionary with parsed attributes (gamma, bootstrap, etc.).
    """
    attrs: dict[str, Any] = {}
    
    if not comment or not comment.startswith('&'):
        return attrs
    
    # Remove leading '&'
    content = comment[1:].strip()
    
    # Parse key=value pairs separated by commas
    parts = content.split(',')
    for part in parts:
        part = part.strip()
        if '=' not in part:
            continue
        
        key, value = part.split('=', 1)
        key = key.strip()
        value = value.strip()
        
        # Try to convert to float
        try:
            attrs[key] = float(value)
        except ValueError:
            attrs[key] = value
    
    return attrs


def from_enewick(enewick_string: str, **kwargs: Any) -> 'DirectedPhyNetwork':
    """
    Parse an eNewick string and create a DirectedPhyNetwork.
    
    This function parses an Extended Newick (eNewick) format string and converts
    it to a DirectedPhyNetwork. It supports:
    - Branch lengths on edges
    - Hybrid nodes (reticulations) using #H markers
    - Gamma and bootstrap values in comments
    - Node labels (quoted and unquoted)
    - Internal node labels
    
    Parameters
    ----------
    enewick_string : str
        The eNewick format string to parse. Must end with ';'.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    DirectedPhyNetwork
        Parsed directed phylogenetic network.
    
    Raises
    ------
    ENewickParseError
        If the eNewick string is malformed or cannot be parsed.
    ValueError
        If the parsed network structure is invalid for DirectedPhyNetwork.
    
    Examples
    --------
    >>> # Simple tree
    >>> net = from_enewick("((A,B),C);")
    >>> net.number_of_nodes()
    4
    >>> net.number_of_edges()
    3
    
    >>> # Tree with branch lengths
    >>> net = from_enewick("((A:0.5,B:0.3):0.1,C:0.2);")
    >>> net.get_branch_length(0, 'A')
    0.5
    
    Notes
    -----
    - Comments starting with '&' are parsed for edge attributes (gamma, bootstrap)
    - Hybrid nodes are identified by #H markers
    - Internal nodes without labels get auto-generated integer IDs
    - Leaves use their label as the node ID
    """
    # Import here to avoid circular dependency
    from .base import DirectedPhyNetwork
    
    # Parse the eNewick string
    parsed = parse_enewick(enewick_string)
    
    # Convert ParsedENewick to DirectedPhyNetwork
    # Build edges list with attributes
    edges: list[dict[str, Any]] = []
    
    # Build a mapping of node IDs to their comment attributes
    node_comments: dict[Any, dict[str, Any]] = {}
    for node in parsed.nodes:
        node_id = node['id']
        if 'comment' in node:
            comment_attrs = _parse_comment_attributes(node['comment'])
            if comment_attrs:
                node_comments[node_id] = comment_attrs
    
    # Process edges
    # Comments on nodes represent attributes for edges TO that node
    for edge in parsed.edges:
        edge_dict: dict[str, Any] = {
            'u': edge['u'],
            'v': edge['v']
        }
        
        # Copy edge attributes (branch_length is already on edge)
        if 'branch_length' in edge:
            edge_dict['branch_length'] = edge['branch_length']
        
        # Check if target node has comment attributes (these are edge attributes for edge TO the node)
        if edge['v'] in node_comments:
            edge_dict.update(node_comments[edge['v']])
        
        # Handle edge key if present
        if 'key' in edge and edge['key'] != 0:
            edge_dict['key'] = edge['key']
        
        edges.append(edge_dict)
    
    # Build nodes list with labels
    nodes: list[tuple[Any, dict[str, Any]]] = []
    for node in parsed.nodes:
        node_id = node['id']
        node_attrs: dict[str, Any] = {}
        
        # Add label if present
        if 'label' in node:
            node_attrs['label'] = node['label']
        
        # Add other node attributes (excluding comment if it was parsed as edge attr)
        for key, value in node.items():
            if key not in ('id', 'label', 'comment'):
                node_attrs[key] = value
        
        # If comment exists but wasn't parsed as edge attributes, keep it
        if 'comment' in node and node['id'] not in node_comments:
            node_attrs['comment'] = node['comment']
        
        nodes.append((node_id, node_attrs))
    
    # Create the network
    return DirectedPhyNetwork(edges=edges, nodes=nodes if nodes else None)

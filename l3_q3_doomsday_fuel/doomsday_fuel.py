# Matrix lib copy from  https://github.com/dmulholland/pymatrix/blob/master/pymatrix.py

import fractions
from fractions import Fraction
import math
import operator
import functools
import textwrap


# Library version.
__version__ = '2.1.1'


# Exports.
__all__ = ['matrix', 'Matrix', 'MatrixError', 'dot', 'cross']


# --------------------------------------------------------------------------
# Library Functions
# --------------------------------------------------------------------------


def matrix(*pargs, **kwargs):
    """ Convenience function for instantiating Matrix objects. """
    if isinstance(pargs[0], int):
        return Matrix.identity(pargs[0])
    elif isinstance(pargs[0], str):
        return Matrix.from_string(*pargs, **kwargs)
    elif isinstance(pargs[0], list):
        return Matrix.from_list(*pargs, **kwargs)
    else:
        raise NotImplementedError


def dot(u, v):
    """ Returns u . v - the scalar product of vectors u and v. """
    return sum(map(operator.mul, u, v))


def cross(u, v):
    """ Returns u x v - the vector product of 3D column vectors u and v. """
    w = Matrix(3, 1)
    w[0][0] = u[1][0] * v[2][0] - u[2][0] * v[1][0]
    w[1][0] = u[2][0] * v[0][0] - u[0][0] * v[2][0]
    w[2][0] = u[0][0] * v[1][0] - u[1][0] * v[0][0]
    return w


# --------------------------------------------------------------------------
# Library Classes
# --------------------------------------------------------------------------


class MatrixError(Exception):
    """ Invalid operation attempted on a Matrix object. """
    pass


class Matrix:

    """ Matrix object supporting basic linear algebra operations. """

    def __init__(self, rows, cols, fill=0):
        """ Initialize a `rows` x `cols` matrix filled with `fill`. """
        self.numrows = rows
        self.numcols = cols
        self.grid = [[fill for i in range(cols)] for j in range(rows)]

    def __str__(self):
        """ Returns a string representation of the matrix. """
        maxlen = max(len(str(e)) for e in self)
        string = '\n'.join(
            ' '.join(str(e).rjust(maxlen) for e in row) for row in self.grid
        )
        return textwrap.dedent(string)

    def __repr__(self):
        """ Returns a string representation of the object. """
        return '<%s %sx%s 0x%x>' % (
            self.__class__.__name__, self.numrows, self.numcols, id(self)
        )

    def __getitem__(self, key):
        """ Enables `self[row][col]` indexing and assignment. """
        return self.grid[key]

    def __contains__(self, item):
        """ Containment: `item in self`. """
        for element in self:
            if element == item:
                return True
        return False

    def __neg__(self):
        """ Negative operator: `- self`. Returns a negated copy. """
        return self.map(lambda element: -element)

    def __pos__(self):
        """ Positive operator: `+ self`. Returns a copy. """
        return self.map(lambda element: element)

    def __eq__(self, other):
        """ Equality: `self == other`. """
        return self.equals(other)

    def __ne__(self, other):
        """ Inequality: `self != other`. """
        return not self.__eq__(other)

    def __add__(self, other):
        """ Addition: `self + other`. """
        if not isinstance(other, Matrix):
            raise MatrixError('cannot add %s to a matrix' % type(other))
        if self.numrows != other.numrows or self.numcols != other.numcols:
            raise MatrixError('cannot add matrices of different sizes')
        m = Matrix(self.numrows, self.numcols)
        for row, col, element in self.elements():
            m[row][col] = element + other[row][col]
        return m

    def __sub__(self, other):
        """ Subtraction: `self - other`. """
        if not isinstance(other, Matrix):
            raise MatrixError('cannot subtract %s from a matrix' % type(other))
        if self.numrows != other.numrows or self.numcols != other.numcols:
            raise MatrixError('cannot subtract matrices of different sizes')
        m = Matrix(self.numrows, self.numcols)
        for row, col, element in self.elements():
            m[row][col] = element - other[row][col]
        return m

    def __mul__(self, other):
        """ Multiplication: `self * other`. """
        if isinstance(other, Matrix):
            if self.numcols != other.numrows:
                raise MatrixError('incompatible sizes for multiplication')
            m = Matrix(self.numrows, other.numcols)
            for row, col, element in m.elements():
                for re, ce in zip(self.row(row), other.col(col)):
                    m[row][col] += re * ce
            return m
        else:
            return self.map(lambda element: element * other)

    def __rmul__(self, other):
        """ Multiplication: `other * self`. Note that this method is never
        called when `other` is a Matrix object - in that case the left
        matrix would handle the multiplication itself via its own __mul__
        method. This method is intended to handle multiplication on the left
        by simple numerical types. """
        return self * other

    def __pow__(self, other):
        """ Exponentiation: `self ** other`. """
        if not isinstance(other, int) or other < 1:
            raise MatrixError('only positive integer powers are supported')
        m = self.copy()
        for i in range(other - 1):
            m = m * self
        return m

    def __iter__(self):
        """ Iteration: `for i in self`. """
        for row in range(self.numrows):
            for col in range(self.numcols):
                yield self[row][col]

    def row(self, n):
        """ Returns an iterator over the specified row. """
        for col in range(self.numcols):
            yield self[n][col]

    def col(self, n):
        """ Returns an iterator over the specified column. """
        for row in range(self.numrows):
            yield self[row][n]

    def rows(self):
        """ Returns a row iterator for each row in the matrix. """
        for row in range(self.numrows):
            yield self.row(row)

    def cols(self):
        """ Returns a column iterator for each column in the matrix. """
        for col in range(self.numcols):
            yield self.col(col)

    def rowvec(self, n):
        """ Returns the specified row as a new row vector. """
        v = Matrix(1, self.numcols)
        for col in range(self.numcols):
            v[0][col] = self[n][col]
        return v

    def colvec(self, n):
        """ Returns the specified column as a new column vector. """
        v = Matrix(self.numrows, 1)
        for row in range(self.numrows):
            v[row][0] = self[row][n]
        return v

    def equals(self, other, delta=None):
        """ Returns true if `self` and `other` are identically-sized matrices
        and their corresponding elements agree to within `delta`. If `delta`
        is omitted, we perform a simple equality check (`==`) on corresponding
        elements instead. """
        if self.numrows != other.numrows or self.numcols != other.numcols:
            return False
        if delta:
            for row, col, element in self.elements():
                if abs(element - other[row][col]) > delta:
                    return False
        else:
            for row, col, element in self.elements():
                if element != other[row][col]:
                    return False
        return True

    def elements(self):
        """ Iterator returning the tuple (row, col, element). """
        for row in range(self.numrows):
            for col in range(self.numcols):
                yield row, col, self[row][col]

    def copy(self):
        """ Returns a copy of the matrix. """
        return self.map(lambda element: element)

    def trans(self):
        """ Returns the transpose of the matrix as a new object. """
        m = Matrix(self.numcols, self.numrows)
        for row, col, element in self.elements():
            m[col][row] = element
        return m

    def det(self):
        """ Returns the determinant of the matrix. """
        if not self.is_square():
            raise MatrixError('non-square matrix does not have determinant')
        ref, _, multiplier = get_row_echelon_form(self)
        ref_det = functools.reduce(
            operator.mul,
            (ref[i][i] for i in range(ref.numrows))
        )
        return ref_det / multiplier

    def minor(self, row, col):
        """ Returns the specified minor. """
        return self.del_row_col(row, col).det()

    def cofactor(self, row, col):
        """ Returns the specified cofactor. """
        return pow(-1, row + col) * self.minor(row, col)

    def cofactors(self):
        """ Returns the matrix of cofactors as a new object. """
        m = Matrix(self.numrows, self.numcols)
        for row, col, element in self.elements():
            m[row][col] = self.cofactor(row, col)
        return m

    def adjoint(self):
        """ Returns the adjoint matrix as a new object. """
        return self.cofactors().trans()

    def inv(self):
        """ Returns the inverse matrix if it exists or raises MatrixError. """
        if not self.is_square():
            raise MatrixError('non-square matrix cannot have an inverse')
        identity = Matrix.identity(self.numrows)
        rref, inverse = get_reduced_row_echelon_form(self, identity)
        if rref != identity:
            raise MatrixError('matrix is non-invertible')
        return inverse

    def del_row_col(self, row_to_delete, col_to_delete):
        """ Returns a new matrix with the specified row & column deleted. """
        return self.del_row(row_to_delete).del_col(col_to_delete)

    def del_row(self, row_to_delete):
        """ Returns a new matrix with the specified row deleted. """
        m = Matrix(self.numrows - 1, self.numcols)
        for row, col, element in self.elements():
            if row < row_to_delete:
                m[row][col] = element
            elif row > row_to_delete:
                m[row - 1][col] = element
        return m

    def del_col(self, col_to_delete):
        """ Returns a new matrix with the specified column deleted. """
        m = Matrix(self.numrows, self.numcols - 1)
        for row, col, element in self.elements():
            if col < col_to_delete:
                m[row][col] = element
            elif col > col_to_delete:
                m[row][col - 1] = element
        return m

    def map(self, func):
        """ Forms a new matrix by mapping `func` to each element. """
        m = Matrix(self.numrows, self.numcols)
        for row, col, element in self.elements():
            m[row][col] = func(element)
        return m

    def rowop_multiply(self, row, m):
        """ In-place row operation. Multiplies the specified row by `m`. """
        for col in range(self.numcols):
            self[row][col] = self[row][col] * m

    def rowop_swap(self, r1, r2):
        """ In-place row operation. Interchanges the two specified rows. """
        for col in range(self.numcols):
            self[r1][col], self[r2][col] = self[r2][col], self[r1][col]

    def rowop_add(self, r1, m, r2):
        """ In-place row operation. Adds `m` times row `r2` to row `r1`. """
        for col in range(self.numcols):
            self[r1][col] = self[r1][col] + m * self[r2][col]

    def ref(self):
        """ Returns the row echelon form of the matrix. """
        return get_row_echelon_form(self)[0]

    def rref(self):
        """ Returns the reduced row echelon form of the matrix. """
        return get_reduced_row_echelon_form(self)[0]

    def len(self):
        """ Vectors only. Returns the length of the vector. """
        return math.sqrt(sum(e ** 2 for e in self))

    def dir(self):
        """ Vectors only. Returns a unit vector in the same direction. """
        return (1 / self.len()) * self

    def is_square(self):
        """ True if the matrix is square. """
        return self.numrows == self.numcols

    def is_invertible(self):
        """ True if the matrix is invertible. """
        try:
            inverse = self.inv()
            return True
        except MatrixError:
            return False

    def rank(self):
        """ Returns the rank of the matrix. """
        rank = 0
        for row in self.ref().rows():
            for element in row:
                if element != 0:
                    rank += 1
                    break
        return rank

    def dot(self, other):
        """ Returns the scalar product: `self` . `other`. """
        return dot(self, other)

    def cross(self, other):
        """ Returns the vector product: `self` x `other`. """
        return cross(self, other)

    @staticmethod
    def from_list(l):
        """ Instantiates a new matrix object from a list of lists. """
        m = Matrix(len(l), len(l[0]))
        for rownum, row in enumerate(l):
            for colnum, element in enumerate(row):
                m[rownum][colnum] = element
        return m

    @staticmethod
    def from_string(s, rowsep=None, colsep=None, parser=fractions.Fraction):
        """ Instantiates a new matrix object from a string. """
        rows = s.strip().split(rowsep) if rowsep else s.strip().splitlines()
        m = Matrix(len(rows), len(rows[0].split(colsep)))
        for rownum, row in enumerate(rows):
            for colnum, element in enumerate(row.split(colsep)):
                m[rownum][colnum] = parser(element)
        return m

    @staticmethod
    def identity(n):
        """ Instantiates a new n x n identity matrix. """
        m = Matrix(n, n)
        for i in range(n):
            m[i][i] = 1
        return m


# --------------------------------------------------------------------------
# Algorithms
# --------------------------------------------------------------------------


# We determine the row echelon form of the matrix using the forward phase of
# the Gauss-Jordan elimination algorithm. If a `mirror` matrix is supplied,
# we apply the same sequence of row operations to it. Note that neither
# matrix is altered in-place; instead copies are returned.
def get_row_echelon_form(matrix, mirror=None):
    matrix = matrix.copy()
    mirror = mirror.copy() if mirror else None
    det_multiplier = 1

    # Start with the top row and work downwards.
    for top_row in range(matrix.numrows):

        # Find the leftmost column that is not all zeros.
        # Note: this step is sensitive to small rounding errors around zero.
        found = False
        for col in range(matrix.numcols):
            for row in range(top_row, matrix.numrows):
                if matrix[row][col] != 0:
                    found = True
                    break
            if found:
                break
        if not found:
            break

        # Get a non-zero entry at the top of this column.
        if matrix[top_row][col] == 0:
            matrix.rowop_swap(top_row, row)
            det_multiplier *= -1
            if mirror:
                mirror.rowop_swap(top_row, row)

        # Make this entry '1'.
        if matrix[top_row][col] != 1:
            multiplier = 1 / matrix[top_row][col]
            matrix.rowop_multiply(top_row, multiplier)
            matrix[top_row][col] = 1 # assign directly in case of rounding errors
            det_multiplier *= multiplier
            if mirror:
                mirror.rowop_multiply(top_row, multiplier)

        # Make all entries below the leading '1' zero.
        for row in range(top_row + 1, matrix.numrows):
            if matrix[row][col] != 0:
                multiplier = -matrix[row][col]
                matrix.rowop_add(row, multiplier, top_row)
                if mirror:
                    mirror.rowop_add(row, multiplier, top_row)

    return matrix, mirror, det_multiplier


# Determine the reduced row echelon form of the matrix using the Gauss-Jordan
# elimination algorithm. If a `mirror` matrix is supplied, the same sequence
# of row operations will be applied to it. Note that neither matrix is
# altered in-place; instead copies are returned.
def get_reduced_row_echelon_form(matrix, mirror=None):

    # Forward phase: determine the row echelon form.
    matrix, mirror, ignore = get_row_echelon_form(matrix, mirror)

    # The backward phase of the algorithm. For each row, starting at the
    # bottom and working up, find the column containing the leading 1 and
    # make all the entries above it zero.
    for last_row in range(matrix.numrows - 1, 0, -1):
        for col in range(matrix.numcols):
            if matrix[last_row][col] == 1:
                for row in range(last_row):
                    if matrix[row][col] != 0:
                        multiplier = -matrix[row][col]
                        matrix.rowop_add(row, multiplier, last_row)
                        if mirror:
                            mirror.rowop_add(row, multiplier, last_row)
                break

    return matrix, mirror


######################################################
#  Solution start from here
####################################################


# gcd, lcm function from https://gist.github.com/endolith/114336
def gcd(*numbers):
    """Return the greatest common divisor of the given integers"""
    from fractions import gcd
    return reduce(gcd, numbers)

# Least common multiple is not in standard libraries? It's in gmpy, but this is simple enough:

def lcm(*numbers):
    """Return lowest common multiple."""
    def lcm(a, b):
        return (a * b) // gcd(a, b)
    return reduce(lcm, numbers, 1)


def preprocess_matrix(origin_matrix):
    size = len(origin_matrix)
    shift_tracker=[]
    prob_matrix = []
    final_count = 0
    is_good_shape = True

    for i in range(size):
        info = {
            'is_final_state': False,
            'new_index': -1,
            'old_index': i
        }

        curr_row = origin_matrix[i]

        row_sum = sum(curr_row)
        if row_sum:
            curr_row = [Fraction(item, row_sum) for item in curr_row]
        else:
            info['is_final_state'] = True
            final_count += 1
            curr_row[i] = 1;
        prob_matrix.append(curr_row)
        shift_tracker.append(info)

    # check if already in good shape
    final_index = size - final_count
    for i in range(final_count):
        if not shift_tracker[final_index]['is_final_state']:
            is_good_shape = False
            break
        final_index += 1
    return prob_matrix, shift_tracker, final_count, is_good_shape


def answer(m):
    # special case
    if len(m) == 1:
        return [1, 1]

    size = len(m)
    prob_mat, shift_tracker, final_count, is_good_shape = preprocess_matrix(m)

    if not is_good_shape:
        temp_mat = [None] * size
        # shift all final state to bottom
        unstable_index = 0
        final_index = size - final_count
        for i in range(size):
            row = prob_mat[i]
            if shift_tracker[i]['is_final_state']:
                temp_mat[final_index] = row
                shift_tracker[i]['new_index'] = final_index
                final_index += 1
            else:
                temp_mat[unstable_index] = row
                shift_tracker[i]['new_index'] = unstable_index
                unstable_index += 1


        # swap columns
        for i in range(final_count):
            # start from bottom right
            idx = size - 1 - i
            if temp_mat[idx][idx] != 1:
                one_index = temp_mat[idx].index(1)
                for j in range(size):
                    temp = temp_mat[j][one_index]
                    temp_mat[j][one_index] = temp_mat[j][idx]
                    temp_mat[j][idx] = temp

        prob_mat = temp_mat

    # shape is good, time to do some matrix operation
    q_mat_size = size - final_count
    final_index = size - final_count
    q_mat = [None] * q_mat_size
    r_mat = [None] * q_mat_size
    # get q, r
    for i in range(q_mat_size):
        q_mat[i] = prob_mat[i][:q_mat_size]
        r_mat[i] = prob_mat[i][final_index:]
    # get i
    i_mat = [0]*q_mat_size
    for i in range(q_mat_size):
        i_mat[i] = [0]*q_mat_size
        i_mat[i][i] = 1

    q_mat = matrix(q_mat)
    i_mat = matrix(i_mat)
    r_mat = matrix(r_mat)

    i_q = i_mat - q_mat
    n_mat = i_q.inv()
    b_mat = n_mat * r_mat

    ns = []
    ds = []
    for item in b_mat.row(0):
        ns.append(item.numerator)
        ds.append(item.denominator)

    denominator = lcm(*ds)

    result = []
    for i in range(len(ns)):
        result.append(ns[i] * denominator / ds[i])
    result.append(denominator)

    return result


# Test
assert (
    answer([
        [0, 2, 1, 0, 0],
        [0, 0, 0, 3, 4],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0]
    ]) == [7, 6, 8, 21]
)

assert (
    answer([
        [0, 1, 0, 0, 0, 1],
        [4, 0, 0, 3, 2, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0]
    ]) == [0, 3, 2, 9, 14]
)

assert (
    answer([
        [1, 2, 3, 0, 0, 0],
        [4, 5, 6, 0, 0, 0],
        [7, 8, 9, 1, 0, 0],
        [0, 0, 0, 0, 1, 2],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0]
    ]) == [1, 2, 3]
)
assert (
    answer([
        [0]
    ]) == [1, 1]
)

assert (
    answer([
        [0, 0, 12, 0, 15, 0, 0, 0, 1, 8],
        [0, 0, 60, 0, 0, 7, 13, 0, 0, 0],
        [0, 15, 0, 8, 7, 0, 0, 1, 9, 0],
        [23, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [37, 35, 0, 0, 0, 0, 3, 21, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ]) == [1, 2, 3, 4, 5, 15]
)

assert (
    answer([
        [0, 7, 0, 17, 0, 1, 0, 5, 0, 2],
        [0, 0, 29, 0, 28, 0, 3, 0, 16, 0],
        [0, 3, 0, 0, 0, 1, 0, 0, 0, 0],
        [48, 0, 3, 0, 0, 0, 17, 0, 0, 0],
        [0, 6, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ]) == [4, 5, 5, 4, 2, 20]
)

assert (
    answer([
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ]) == [1, 1, 1, 1, 1, 5]
)

assert (
    answer([
        [1, 1, 1, 0, 1, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 1, 1, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 1, 1, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 1, 0, 1, 1, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 1, 0, 1, 0, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ]) == [2, 1, 1, 1, 1, 6]
)

assert (
    answer([
        [0, 86, 61, 189, 0, 18, 12, 33, 66, 39],
        [0, 0, 2, 0, 0, 1, 0, 0, 0, 0],
        [15, 187, 0, 0, 18, 23, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ]) == [6, 44, 4, 11, 22, 13, 100]
)

assert (
    answer([
        [0, 0, 0, 0, 3, 5, 0, 0, 0, 2],
        [0, 0, 4, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 4, 4, 0, 0, 0, 1, 1],
        [13, 0, 0, 0, 0, 0, 2, 0, 0, 0],
        [0, 1, 8, 7, 0, 0, 0, 1, 3, 0],
        [1, 7, 0, 0, 0, 0, 0, 2, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ]) == [1, 1, 1, 2, 5]
)

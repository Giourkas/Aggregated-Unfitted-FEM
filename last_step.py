"""The last step is to interpolate vertices on differnet mesh types, e.g. 2D squares, 2D simplices, 3D simplices, 3D squares.
The 3D case is not tested.
"""
# ngsolve stuff
from ngsolve import *
# Basic geometry features (for the background mesh)
from netgen.geom2d import SplineGeometry
# Visualization
from ngsolve.internal import *
# Basic xfem functionality
from xfem import *
from enum import Enum
import numpy as np

class MeshType(Enum):
	squares_2D = 0
	simplices_2D = 1
	squares_3D = 2
	simplices_3D = 3

class VertexCoefficient:
	"""Class that contains information about which coefficients are
	associated to which vertex when interpolating a certain vertex.

	Attributes:
	
		vertex_id (int):
			the index of the Vertex

		coefficient (float):
			the associated coefficent

	"""
	def __init__(self, vertex_id, coefficient):
		"""Constrtucts a VertexCoefficient"""
		self.vertex_id = vertex_id
		self.coefficient = coefficient
	def __str__(self):
		return "vertex(id:" + str(self.vertex_id) + ",cof:" + str(self.coefficient)+ ")"

class LastStep:
	"""Class to interpolate vertices
	"""
	def __init__(self, mesh, cut_info):
		"""Constructor

		Args:
			mesh (Mesh):
				a netgen 2d mesh object for which we want to compute
			cut_info (CutInfo):
				the cut info to which to compute the clusters.

		"""
		self.mesh = mesh
		self.cut_info = cut_info
		self.mesh_type = MeshType.simplices_2D

		is_3D = (mesh.dim == 3)
		num_corners = len(self.mesh[NodeId(FACE,0)].vertices)
		if is_3D:
			if num_corners == 4:
				self.mesh_type = MeshType.simplices_3D
			if num_corners == 8:
				self.mesh_type = MeshType.squares_3D
		else:
			if num_corners == 3:
				self.mesh_type = MeshType.simplices_2D
			if num_corners == 4:
				self.mesh_type = MeshType.squares_2D

	def element_corners(self, e_id):
		if self.mesh.dim == 3:
			return self.mesh[NodeId(CELL, e_id)].vertices
		else:
			return self.mesh[NodeId(FACE, e_id)].vertices

	def get_coefficients(self, vertex_id, face_id):
		"""
		Calculates the coeffients to interpolate the vertex using an H1
		function of order 1 on the given Element.
		Assumes a 2D space with square Elements.

		Args:

			vertex_id (int):
				The index of the vertex that will be interpolated

			face_id (int): the index of the element where the H1 function is on

		Returns:
			list of VertexCoefficient. 

		"""
		if self.mesh_type == MeshType.simplices_3D:
			return self.do_simplices_3D(vertex_id, face_id)
		if self.mesh_type == MeshType.squares_3D:
			return self.do_squares_3D(vertex_id, face_id)
		if self.mesh_type == MeshType.simplices_2D:
			return self.do_simplices_2D(vertex_id, face_id)
		if self.mesh_type == MeshType.squares_2D:
			return self.do_squares_2D(vertex_id, face_id)

	def do_squares_2D(self, vertex_id, face_id):
		"""
		Calculates the coeffients to interpolate the vertex using an H1
		function of order 1 on the given Element.
		Assumes a 2D space with square Elements.

		Args:
			vertex_id (int):
				The index of the vertex that will be interpolated

			face_id (int):
				The index of the element where the H1 function is on

		Returns:
			list of VertexCoefficient.

		"""
		return self.do_squares(vertex_id, face_id, 2)

	def do_simplices_2D(self, vertex_id, face_id):
		"""
		Claulates the coeffients to interpolate the vertex using a H1
		function of order 1 on the given Element.
		Assumes a 2D space with triangular Elements.

		Args:
			vertex_id (int):
				The index of the vertex that will be interpolated

			face_id (int):
				The index of the element where the H1 function is on

		Returns:
			list of VertexCoefficient. 

		"""
		def one_corner(p_dest, point1, point2, point3):
			mat = np.array([[1, point1[0], point1[1]],
			           [1, point2[0], point2[1]], \
		   	           [1, point3[0], point3[1]], \
			])
			b = np.array([1,0, 0])
			res = np.linalg.solve(mat, b)
			return res[0] + p_dest[0] * res[1] + p_dest[1] * res[2]
		vertex_loc = self.mesh[NodeId(VERTEX, vertex_id)].point
		corners = [self.mesh[v] for v in self.element_corners(face_id)]
		return [ \
			VertexCoefficient(corners[0].nr, one_corner(vertex_loc, corners[0].point, corners[1].point, corners[2].point)), \
			VertexCoefficient(corners[1].nr, one_corner(vertex_loc, corners[1].point, corners[0].point, corners[2].point)), \
			VertexCoefficient(corners[2].nr, one_corner(vertex_loc, corners[2].point, corners[1].point, corners[0].point)), \
		]

	def do_squares(self, vertex_id, face_id, ndims):
		n_corners = 2**ndims
		n_corners_half = 2**(ndims - 1)
		corners = [self.mesh[v] for v in self.element_corners(face_id)]
		vertex_loc = self.mesh[NodeId(VERTEX, vertex_id)].point

		#validate input, in particular the algorithm assumes
		#the vertices to be listed in a certain order
		def is_opposite_corner(v1, v2):
			for a1, a2 in zip(v1, v2):
				if a1 == a2:
					return False
			return True
		assert len(corners) == n_corners
		for i in range(n_corners_half):
			assert is_opposite_corner(corners[i].point, corners[i + n_corners_half].point)

		res = []
		for i in range(n_corners):
			i_opp = (i + n_corners_half) % n_corners
			point = corners[i].point
			point_opp = corners[i_opp].point
			coeff = 1
			for dim in range(ndims):
				coeff *= ((vertex_loc[dim] - point_opp[dim])/(point[dim] - point_opp[dim]))
			res.append(VertexCoefficient(corners[i].nr, coeff))
		return res

	def do_squares_3D(self, vertex_id, face_id):
		return self.do_squares(vertex_id, face_id, 3)

	
	def do_simplices_3D(self, vertex_id, face_id):
		def one_corner(p_dest, point1, point2, point3, point4):
			# f(x) =
			mat = np.array([[1, point1[0], point1[1], point1[2]],
			           [1, point2[0], point2[1], point2[2]], \
		   	           [1, point3[0], point3[1], point3[2]], \
		   	           [1, point4[0], point4[1], point4[2]], \
			])
			b = np.array([1, 0, 0, 0])
			res = np.linalg.solve(mat, b)
			return res[0] + p_dest[0] * res[1] + p_dest[1] * res[2] + p_dest[2] * res[3]
		vertex_loc = self.mesh[NodeId(VERTEX, vertex_id)].point
		corners = [self.mesh[v] for v in self.element_corners(face_id)]
		return [ \
			VertexCoefficient(corners[0].nr, one_corner(vertex_loc, corners[0].point, corners[3].point, corners[1].point, corners[2].point)), \
			VertexCoefficient(corners[1].nr, one_corner(vertex_loc, corners[1].point, corners[3].point, corners[0].point, corners[2].point)), \
			VertexCoefficient(corners[2].nr, one_corner(vertex_loc, corners[2].point, corners[3].point, corners[1].point, corners[0].point)), \
			VertexCoefficient(corners[2].nr, one_corner(vertex_loc, corners[3].point, corners[1].point, corners[2].point, corners[0].point)), \
		]

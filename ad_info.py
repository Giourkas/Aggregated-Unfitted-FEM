"""The ad_info module below gives a new class with information about adjacent faces.	Specifically it gives info about:

a) The connection of vertices and faces.

b) A list of adjacent faces.

c) Adjacent faces of a vertex.

d) Two test functions to see in the mesh any adjacent faces or adjacent faces of a vertex.
"""

# ngsolve stuff
from ngsolve import *
# basic geometry features (for the background mesh)
from netgen.geom2d import SplineGeometry
# visualization stuff
from ngsolve.internal import *
# basic xfem functionality
from xfem import *

class AdInfo:
	'''
		Information about which faces are adjacent.
	'''
	def __init__(self, mesh):
		"""
		Args:
			mesh (Mesh from netgen):
		"""
		#Proerties of the AdInfo structure
		self.mesh = mesh
	
		nelements = self.mesh.nnodes(self.get_element_name())
		self.element_to_face = []
		self.face_to_element = [0 for x in range(nelements)]
		if self.mesh.dim == 2:
			self.element_to_face = [e.faces[0].nr for e in self.mesh.Elements(VOL)]
		else:
			self.element_to_face = [e.nr for e in self.mesh.Elements(VOL)]
		for e_id, f_id in enumerate(self.element_to_face):
			self.face_to_element[f_id] = e_id

		self.f2f = []
		self.v2f = []
		
		#Initilaize  self.f2f and self.v2f
		self.build_f2f()
		self.build_v2f()


	def get_facet_name(self):
		return EDGE if self.mesh.dim == 2 else FACE;

	def get_facets(self, e_id):
		element = self.mesh[NodeId(self.get_element_name(), e_id)]
		return element.edges if self.mesh.dim == 2 else element.faces;

	def get_element_name(self):
		return FACE if self.mesh.dim == 2 else CELL;

	def build_f2f(self):
		"""Gathers information about connection between edges to faces first,
		and faces to faces afterwards in the mesh. A helper function for the
		constructor, don't call manually.
		"""
		mesh = self.mesh
		nelements = self.mesh.nnodes(self.get_element_name())
		nfacets = self.mesh.nnodes(self.get_facet_name())
		#Assumes a 2d mesh.
		#A map that maps every edge to its two adjacent faces.
		#Used to build face_to_face.
		edge_to_face = [[] for x in range(nfacets)]
		#A map that maps every face to its (at most 3 for 2d simplices) adjacent faces.
		self.f2f = [[] for x in range(nelements)]
		#Build edge_to_face
		#For f in mesh.faces:
		for e_id in range(nelements):
			#print("face nr", f.nr)
			for facet in self.get_facets(e_id):
				edge_to_face[facet.nr].append(e_id)
		#Build self.f2f
		for e in edge_to_face:
			#Ignore edges that are not adjacent to two faces (border edges)
			if len(e) > 1:
				self.f2f[e[0]].append(e[1])
				self.f2f[e[1]].append(e[0])

	def build_v2f(self):
		"""Gathers information about connection between verticies to faces in the
		mesh. A helper function for the constructor, don't call manually.
		"""
		mesh = self.mesh
		nelements = self.mesh.nnodes(self.get_element_name())
		self.v2f = [[] for x in range(self.mesh.nnodes(VERTEX))]
		for e_id in range(nelements):
			element = self.mesh[NodeId(self.get_element_name(), e_id)]
			for vertex in element.vertices:
				self.v2f[vertex.nr].append(e_id)

	def get_ajacent_faces(self, face_index):
		"""Function that gives list of adjacent faces of a face.
		Args:
		face_index (int): Index of the face
			
		Returns: Adjacent faces.
		"""
		return self.f2f[face_index]

	def get_vertex_ajacent_faces(self, vertex_index):
		"""Function for the adjacent faces of a vertex.
		Args:
		vertex_index (int): Index of the vertex 
			
		Returns: Adjacent faces of the vertex.
		"""
		return self.v2f[vertex_index]

	def test_f2f(self, face_index):
		"""Test function to visualize the result of the get_ajacent_faces function.
		Args:
		face_index (int): Index of a face

		"""
		l2 = L2(self.mesh)
		gf_e = GridFunction(l2)
		for f_id in self.get_ajacent_faces(face_index):
			gf_e.vec[self.face_to_element[f_id]] = True
		Draw(gf_e, self.mesh,"test_f2f_" + str(face_index))

	def test_v2f(self, vertex_index):
		"""Test function to visualize the result of the get_vertex_ajacent_faces.
		Args:
		vertex_index (int): Index of a vertex

		"""
		l2 = L2(self.mesh)
		gf_e = GridFunction(l2)
		for f_id in self.get_vertex_ajacent_faces(vertex_index):
			gf_e.vec[self.face_to_element[f_id]] = True
		Draw(gf_e, self.mesh,"test_v2f_" + str(vertex_index))

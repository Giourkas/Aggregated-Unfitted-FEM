"""This is the basic module of the algorithm to build the clusters. Firstly, it defines a new class to generate the clusters and another function to give the distance between the cells. The idea is to build a cluster between an inner cell and neighbour cut cells having the smallest distance. There are also two functions defined for the visualization of clusters created.
"""
# ngsolve stuff
from ngsolve import *
# Basic geometry features (for the background mesh)
from netgen.geom2d import SplineGeometry
# Visualization
from ngsolve.internal import *
# Basic xfem functionality
from xfem import *
import itertools
import numpy as np
import time

def copy_bitarray(a):
	"""Creates a copy of a BitArray structure

	Args:
		a (BitArray):

	Returns:
		BitArray, A copy of a

	"""
	copy = BitArray(len(a))
	for i in range(len(a)):
		copy[i] = a[i]
	return copy

class CusterInfo:
	"""Class to generate the clusters when doing unfitted FEM
	"""
	def __init__(self, ad_info, cut_info):
		"""Constructor

		Arguments:
			ad_info (AdInfo):
				Usually created via AdInfo(mesh), where mesh is cut_info.Mesh()

			cut_info (CutInfo):
				The cut info to which to compute the clusters.
		"""
		self.mesh = cut_info.Mesh()
		self.ad_info = ad_info
		self.cut_info = cut_info
		self.cluster_map = {}
		self.distance_map = {}
		self.build_clusters()
		self.ci_pos = self.GetFacesOfType(POS)
		self.ci_neg = self.GetFacesOfType(NEG)
		self.ci_if = self.GetFacesOfType(IF)

	def element_corners(self, e_id):
		"""returns a list of all verticies adjacent to an element
		"""
		if self.mesh.dim == 3:
			return self.mesh[NodeId(CELL, e_id)].vertices
		else:
			return self.mesh[NodeId(FACE, e_id)].vertices

	def distance_cell_to_cell(self, f_id1, f_id2):
		"""Helper method to compute the distance of two elements, just
		computes the distance between the middlepoints of the two elements.

		Args:
			f_id1 (int):
				The index of the first Element

			f_id2 (int):
				The index of the second Element

		Returns:
			float. The distance of the two Elements
		"""
		corners1 = [np.asarray(self.mesh[v].point) for v in self.element_corners(f_id1)]
		corners2 = [np.asarray(self.mesh[v].point) for v in self.element_corners(f_id2)]
		return np.linalg.norm(sum(corners1)/len(corners1) - sum(corners2)/len(corners2))

	def GetFacesOfType(self, t):
		"""wrapper for cut_info.GetElementsOfType that uses the same
		ordering as used by mesh[NodeId(FACE, id)]
		"""
		el = self.cut_info.GetElementsOfType(t)
		res = BitArray(len(el))
		for i in range(len(el)):
			res[self.ad_info.element_to_face[i]] =  el[i]
		return res

	def build_clusters(self):
		"""Helper method to do the actual work, only called by the constructor
		of this class
		"""
		self.cluster_map = {}
		self.distance_map = {}
		ci_neg = self.GetFacesOfType(NEG)
		ci_if = self.GetFacesOfType(IF)
		#Important: ci_neg, ci_if are direct references to the ci object, so we may not change them
		#So we make a copy here.
		#Use a bitset instead of a list of indices so that we change it while iterating without problems.
		faces_left = copy_bitarray(ci_if)

		#The first iteration, add all split faces to adjacent inner faces
		for f_id, is_in in enumerate(faces_left):
			if is_in:
				#Assign border faces of adjacent inner faces
				candidates = [ad_f_id for ad_f_id in self.ad_info.get_ajacent_faces(f_id) if ci_neg[ad_f_id]]
				if False:
					pass
				else:
					for ad_f_id in self.ad_info.get_ajacent_faces(f_id):
						if ci_neg[ad_f_id]:
							self.cluster_map[f_id] = ad_f_id
							faces_left[f_id] = False
							self.distance_map[f_id] = 1

		for i in itertools.count(1):
			# Once all border faces are assigned we are done.
			left_todo = faces_left.NumSet()
			if left_todo == 0:
				break
			for f_id, is_in in enumerate(faces_left):
				if is_in:
					# Assign the face to the same inner face as any adjacent border cell.
					candidates = [ad_f_id for ad_f_id in self.ad_info.get_ajacent_faces(f_id) if ad_f_id in self.cluster_map and self.distance_map[ad_f_id] == i]
					if len(candidates) > 0:
						nearest_f_id = min(candidates, key=lambda ad_face_id: self.distance_cell_to_cell(ad_face_id, f_id))
						faces_left[f_id] = False
						self.cluster_map[f_id] = self.cluster_map[nearest_f_id]
						self.distance_map[f_id] = i + 1
					if False:
						for ad_f_id in self.ad_info.get_ajacent_faces(f_id):
							if ad_f_id in self.cluster_map and self.distance_map[ad_f_id] == i:
								faces_left[f_id] = False
								self.cluster_map[f_id] = self.cluster_map[ad_f_id]
								self.distance_map[f_id] = i + 1
			# Assert that we make progress. This can only fail if there is a cut
			# face that is not connected to any inner face
			# assert(left_todo != faces_left.NumSet())

	def get_cluster(self, face_id):
		"""Returns the id of the cluster the Element belongs to, that id
		is the same as the Element id of the single uncut element in that cluster.

		Args:
			face_id (int):
				The index of the element

		Returns:
			int. The index of the cluster.
		"""
		assert (not self.ci_pos[face_id]), "outer faces have no cluster"
		if self.ci_neg[face_id]:
			return face_id
		else:
			return self.cluster_map[face_id]

	def draw_clusters_inone(self):
		"""Tries to draw all clusters in one image. Due to the limited
		number of colors, this will not always work very well.
		"""
		def shuffle(n):
			return n % 7+ 1
		gf_enumbers = GridFunction(L2(self.mesh))

		for i in range(len(gf_enumbers.vec)):
			f_id = self.ad_info.element_to_face[i]
			if self.ci_pos[f_id]:
				gf_enumbers.vec[i] = 0
			else:
				gf_enumbers.vec[i] = shuffle(self.get_cluster(f_id))
		Draw(gf_enumbers, self.mesh, "clusters")

	def draw_cluster(self, f_id):
		"""Draws the cluster with the given id, that is the id of the
		unique uncut element in that cluster

		Args:
			face_id (int):
				The index of the element
		"""
		gf_enumbers = GridFunction(L2(self.mesh))
		count = 0
		for i in range(len(gf_enumbers.vec)):
			f2_id = self.ad_info.element_to_face[i]
			if f_id == f2_id:
				count += 1
				gf_enumbers.vec[i] = 1
			elif self.ci_pos[f2_id]:
				gf_enumbers.vec[i] = 0
			elif self.get_cluster(f2_id) == f_id:
				gf_enumbers.vec[i] = 2
				count += 1
			elif self.ci_neg[f2_id]:
				gf_enumbers.vec[i] = 3
			else:
				gf_enumbers.vec[i] = 4
		if count > 1:
			Draw(gf_enumbers, self.mesh, "cluster" + str(f_id) + "(" + str(count) + ")")

	def draw_clusters(self):
		"""Draws all calculated clusters, each in a diffent image
		named '"cluster<custer_id>(<custer_size>)'
		"""
		for i in range(self.mesh.nnodes(ELEMENT)):
			if self.ci_neg[i]:
				self.draw_cluster(i)

	def get_posible_interpolation_faces(self, vertex_id):
		"""Returns a list of all Elements that could be used to interpolate the
		DOF associated to the vertex with the given id.

		Args:
			vertex_id (int):
				The index of the vertex
		
		Returns:
			list of int. List of indices of Elements.
		"""
		ad_faces = self.ad_info.get_vertex_ajacent_faces(vertex_id)
		assert not any(self.ci_neg[i] for i in ad_faces), "No interpolation needed"
		return [self.get_cluster(face_id) for face_id in ad_faces if self.ci_if[face_id]]

	def get_interpolation_face(self, vertex_id):
		"""Returns the Elements that should be used to interpolate the given vertexs

		Args:
			vertex_id (int):
				The index of the vertex
		
		Returns:
			int. Index of the element.
		"""
		def distance_cell_to_vertex(face_id, vertex_id):
			face_corners = [np.asarray(self.mesh[v].point) for v in self.mesh[NodeId(self.ad_info.get_element_name(), face_id)].vertices]
			face_pos = sum(face_corners)/len(face_corners)
			vertex_pos = np.asarray(self.mesh[NodeId(VERTEX, vertex_id)].point)
			return np.linalg.norm(face_pos-vertex_pos)
		return min(self.get_posible_interpolation_faces(vertex_id), key=lambda face_id: distance_cell_to_vertex(face_id, vertex_id))



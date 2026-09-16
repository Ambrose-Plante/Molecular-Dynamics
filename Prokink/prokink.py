import numpy as np
import math
import mdtraj as md


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def project_onto_plane(vec, plane_normal):

    # vector n: n is orthogonal vector to Plane P 
    n = plane_normal        

    # finding norm of the vector n  
    n_norm = np.sqrt(sum(n**2))     

    # project vector onto the orthogonal vector n 
    # find dot product using np.dot() 
    proj_of_u_on_n = (np.dot(vec, n)/n_norm**2)*n 

    # return projection of u on Plane P 
    return vec - proj_of_u_on_n


def prokink(proline_resid, pre_pro_resid1, pre_pro_resid2, post_pro_resid1, post_pro_resid2, pdb=False, psf=False, dcd=False, checker=True):

    if checker:
        if (post_pro_resid2+post_pro_resid1) <= (pre_pro_resid1+pre_pro_resid2):
            raise Exception("The sum of the post_pro_resid's are less than the sum of the pre_pro_resid's. Check your definitions of pre and post helix :)")
    
        if pre_pro_resid2 <= pre_pro_resid1:
            raise Exception("pre_pro_resid2 must be greater than pre_pro_resid1 :)")
            
        if post_pro_resid2 <= post_pro_resid1:
            raise Exception("post_pro_resid2 must be greater than post_pro_resid1 :)")
            
    if pdb==False:
        traj = md.load(dcd, top=psf)
    else:
        traj = md.load_pdb(pdb)

    # create selection query based on this resid
    selection_query = 'protein and resSeq ' + str(proline_resid) + ' and name CA'

    # get the index of the alpha carbon of this proline
    Pro_CA_index = traj.topology.select(selection_query)

    # get the xyz coordinates of that atom
    pro_CA_xyz = traj.xyz[:, Pro_CA_index[0],:]
    pro_CA_xyz = np.expand_dims(pro_CA_xyz, axis=1)

    # center on that atom in all frames
    # this is important for the defintion of the pre_helix axis
    centered_traj = (traj.xyz - pro_CA_xyz)

    # create selection query for pre-proline helix
    selection_query = 'protein and backbone and resSeq ' + str(pre_pro_resid1) + ' to ' + str(pre_pro_resid2)

    # get the indices of the pre_helix
    pre_helix_indices = traj.topology.select(selection_query)

    # create selection query for post-helix
    selection_query = 'protein and backbone and resSeq ' + str(post_pro_resid1) + ' to ' + str(post_pro_resid2)

    # get the indices of the post-helix
    post_helix_indices = traj.topology.select(selection_query)

    # create selection query for pre-proline helix CA1
    selection_query = 'protein and resSeq ' + str(pre_pro_resid1) + ' and name CA'
    # get the indices
    pre_res1 = traj.topology.select(selection_query)
    # create selection query for pre-proline helix CA2
    selection_query = 'protein and resSeq ' + str(pre_pro_resid2) + ' and name CA'
    # get the indices
    pre_res2 = traj.topology.select(selection_query)

    # create selection query for post-proline helix CA1
    selection_query = 'protein and resSeq ' + str(post_pro_resid1) + ' and name CA'
    # get the indices
    post_res1 = traj.topology.select(selection_query)
    # create selection query for pre-proline helix CA2
    selection_query = 'protein and resSeq ' + str(post_pro_resid2) + ' and name CA'
    # get the indices
    post_res2 = traj.topology.select(selection_query)

    bend_angles = np.zeros(centered_traj.shape[0])
    wobble_angles = np.zeros(centered_traj.shape[0])

    for i,frame in enumerate(centered_traj):

        # BEND ANGLE
        ##############################################
        data = frame[post_helix_indices, :]

        # Calculate the mean of the points, i.e. the 'center' of the cloud
        datamean = data.mean(axis=0)

        # Do an SVD on the mean-centered data.
        uu, dd, vv = np.linalg.svd(data - datamean)

        # The eigenvector coorsponding to the biggest eigenvalue (always the first) is the best fit slope
        post_helix_vec = vv[0]

        ## It is possible that the helix axiz is negated, we must check:
        ## It should make an acute angle with a vector from the first to the last CA in that helix
        # get coordinates of
        res1 = frame[post_res1, :]
        res2 = frame[post_res2, :]

        #if angle is obtuse, flip the vector
        check_angle = math.degrees(angle_between(post_helix_vec, (res2-res1)[0]))
        if check_angle>90:
            post_helix_vec = -post_helix_vec

        # Repeat this procedure for pre_helix
        data = frame[pre_helix_indices, :]
        datamean = data.mean(axis=0)
        uu, dd, vv = np.linalg.svd(data - datamean)
        pre_helix_vec = vv[0]
        # note res definitions are flipped for pre-helix
        res1 = frame[pre_res2, :]
        res2 = frame[pre_res1, :]
        check_angle = math.degrees(angle_between(pre_helix_vec, (res2-res1)[0]))
        if check_angle>90:
            pre_helix_vec = -pre_helix_vec

        # get the angle
        angle = math.degrees(angle_between(pre_helix_vec, post_helix_vec))
        bend_angles[i] = angle
        # FINISH BEND ANGLE
        ########################################################

        # WOBBLE ANGLE
        ##############################################

        # center of pre_helix
        pre_helix_center = frame[pre_helix_indices, :].mean(axis=0)

        # re-center frame to center of pre-helix
        frame = frame - pre_helix_center

        # proline coordinates
        pro_CA_xyz = frame[Pro_CA_index, :].mean(axis=0)

        # vector from pre_helix center to proline alpha carbon (same as coordinates)
        pro_helix_vec = pro_CA_xyz

        # normal vector defining the plane bisecting the pre_helix
        # same as vector defining axis of pre_helix
        ## pre_helix_vec

        # calculate rotation matrix from pre_helix axis to x axis
        mat = rotation_matrix_from_vectors(pre_helix_vec, [1,0,0])

        # rotate entire frame to align pre_helix axis with x zxis
        frame = mat.dot(frame.T).T

        ##### re-calculate all vectors after a rotation of system
        data = frame[post_helix_indices, :]
        datamean = data.mean(axis=0)
        uu, dd, vv = np.linalg.svd(data - datamean)
        post_helix_vec = vv[0]
        res1 = frame[post_res1, :]
        res2 = frame[post_res2, :]
        check_angle = math.degrees(angle_between(post_helix_vec, (res2-res1)[0]))
        if check_angle>90:
            post_helix_vec = -post_helix_vec

        data = frame[pre_helix_indices, :]
        datamean = data.mean(axis=0)
        uu, dd, vv = np.linalg.svd(data - datamean)
        pre_helix_vec = vv[0]
        res1 = frame[pre_res2, :]
        res2 = frame[pre_res1, :]
        check_angle = math.degrees(angle_between(pre_helix_vec, (res2-res1)[0]))
        if check_angle>90:
            pre_helix_vec = -pre_helix_vec

        pro_CA_xyz = frame[Pro_CA_index, :].mean(axis=0)
        pro_helix_vec = pro_CA_xyz
        ###################################

        # component of proline CA vector parellel to YZ plane
        temp = project_onto_plane(pro_helix_vec, pre_helix_vec)

        ###########################
        # Rotate proline CA to Y-axis
        ########################

        # calculate rotation matrix from component of proline CA vector parallel to YZ plane to Y axis
        mat = rotation_matrix_from_vectors(temp, [0,1,0])

        # rotate entire frame to align proline CA with Y zxis
        frame = mat.dot(frame.T).T

        ##### re-calculate all vectors after a rotation of system
        data = frame[post_helix_indices, :]
        datamean = data.mean(axis=0)
        uu, dd, vv = np.linalg.svd(data - datamean)
        post_helix_vec = vv[0]
        res1 = frame[post_res1, :]
        res2 = frame[post_res2, :]
        check_angle = math.degrees(angle_between(post_helix_vec, (res2-res1)[0]))
        if check_angle>90:
            post_helix_vec = -post_helix_vec

        data = frame[pre_helix_indices, :]
        datamean = data.mean(axis=0)
        uu, dd, vv = np.linalg.svd(data - datamean)
        pre_helix_vec = vv[0]
        res1 = frame[pre_res2, :]
        res2 = frame[pre_res1, :]
        check_angle = math.degrees(angle_between(pre_helix_vec, (res2-res1)[0]))
        if check_angle>90:
            pre_helix_vec = -pre_helix_vec

        pro_CA_xyz = frame[Pro_CA_index, :].mean(axis=0)
        pro_helix_vec = pro_CA_xyz
        ###################################

        # component of proline CA vector parellel to YZ plane
        temp = project_onto_plane(pro_helix_vec, pre_helix_vec)

        # component of proline CA vector normal to YZ plane
        temp2 = pro_helix_vec - temp

        # move proline CA to YZ plane
        frame = frame - temp2

        # projection of proline vector onto the plane
        wobble_ref = project_onto_plane(pro_helix_vec, pre_helix_vec)

        # projection of post_helix vector onto the plane
        wobble_post = project_onto_plane(post_helix_vec, pre_helix_vec)

        wobble_angle = math.degrees(angle_between(wobble_ref, wobble_post))

        if wobble_post[2]<0:
            wobble_angle = -wobble_angle

        wobble_angles[i]=wobble_angle

        # FINISH WOBBLE ANGLE
        ########################################################
        
    return bend_angles, wobble_angles


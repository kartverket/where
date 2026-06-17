""" Tests for the data.position module when combining quasar and satellite observations

When combining quasar and satellite observations it is important that the other attribute is of type DirectionArray
and that the other_2 attribute is of type PositionArray or PosVelArray. Changing the order will not work. 
"""


def test_elevation(site_pos, sat_pos, src_dir):
    e_1_sat = site_pos.elevation_to(sat_pos)[0]
    e_1_src = site_pos.elevation_to(src_dir)[1]
    
    #site_pos.other = src_dir
    #site_pos.other_2 = sat_pos
    
    e_2_sat = site_pos.elevation[0]
    e_2_src = site_pos.elevation[1]
    
    assert e_1_sat == e_2_sat
    assert e_1_src == e_2_src  

def test_azimuth(site_pos, sat_pos, src_dir):
    a_1_sat = site_pos.azimuth_to(sat_pos)[0]
    a_1_src = site_pos.azimuth_to(src_dir)[1]
    
    #site_pos.other = src_dir
    #site_pos.other_2 = sat_pos
    
    a_2_sat = site_pos.azimuth[0]
    a_2_src = site_pos.azimuth[1]
    
    assert a_1_sat == a_2_sat
    assert a_1_src == a_2_src

def test_zenith_distance(site_pos, sat_pos, src_dir):
    z_1_sat = site_pos.zenith_distance_to(sat_pos)[0]
    z_1_src = site_pos.zenith_distance_to(src_dir)[1]
    
    #site_pos.other = src_dir
    #site_pos.other_2 = sat_pos
    
    z_2_sat = site_pos.zenith_distance[0]
    z_2_src = site_pos.zenith_distance[1]
    
    assert z_1_sat == z_2_sat
    assert z_1_src == z_2_src

def test_distance(site_pos, sat_pos, src_dir):
    d_1_sat = site_pos.distance_to(sat_pos)[0]
    d_1_src = site_pos.distance_to(src_dir)[1]

    #site_pos.other = src_dir
    #site_pos.other_2 = sat_pos

    d_2_sat = site_pos.distance[0]
    d_2_src = site_pos.distance[1]

    assert d_1_sat == d_2_sat
    assert d_1_src == d_2_src

def test_vector(site_pos, sat_pos, src_dir):
    v_1_sat = site_pos.vector_to(sat_pos)[0]
    v_1_src = site_pos.vector_to(src_dir)[1]

    #site_pos.other = src_dir
    #site_pos.other_2 = sat_pos

    v_2_sat = site_pos.vector[0]
    v_2_src = site_pos.vector[1]

    assert (site_pos.gcrs.vector[0] == (sat_pos.gcrs.pos.val[0] - site_pos.gcrs.pos.val[0])).all()
    assert (site_pos.gcrs.vector[1] == src_dir.gcrs[1].val).all()
    assert (v_1_sat == v_2_sat).all()
    assert (v_1_src == v_2_src).all()


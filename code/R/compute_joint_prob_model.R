# this script computes the joint probability of preferring a specific density-dependent model or geosse by turning on/off
# features, based on two independent networks analyses on anolis and rescaled viburnum tree 

# Network 1: anolis
p_w = 9.75721538066864E-01
p_e = 9.71086919307709E-01
p_d = 6.70620977878571E-01

q_w = 1 - p_w
q_e = 1 - p_e
q_d = 1 - p_d

m_0 = q_w * q_e * q_d
m_w = p_w * q_e * q_d
m_e = q_w * p_e * q_d
m_d = q_w * q_e * p_d
m_we = p_w * p_e * q_d
m_wd = p_w * q_e * p_d
m_ed = q_w * p_e * p_d
m_wed = p_w * p_e * p_d

m_all = c( m_0, m_w, m_e, m_d, m_we, m_wd, m_ed, m_wed )
names(m_all) = c( "0", "w", "e", "d", "we", "wd", "ed", "wed")
print(m_all)

# Network 2: anolis
p_w = 9.96811807155609E-01
p_e = 9.35789406299591E-01
p_d = 8.17522883415222E-01

q_w = 1 - p_w
q_e = 1 - p_e
q_d = 1 - p_d

m_0 = q_w * q_e * q_d
m_w = p_w * q_e * q_d
m_e = q_w * p_e * q_d
m_d = q_w * q_e * p_d
m_we = p_w * p_e * q_d
m_wd = p_w * q_e * p_d
m_ed = q_w * p_e * p_d
m_wed = p_w * p_e * p_d

m_all_2 = c( m_0, m_w, m_e, m_d, m_we, m_wd, m_ed, m_wed )
names(m_all_2) = c( "0", "w", "e", "d", "we", "wd", "ed", "wed")
print(m_all_2)

# Network 1: rescaled viburnum
p_w = 9.4416332244873E-01
p_e = 1.69942229986191E-01
p_d = 8.87794137001038E-01

q_w = 1 - p_w
q_e = 1 - p_e
q_d = 1 - p_d

m_0 = q_w * q_e * q_d
m_w = p_w * q_e * q_d
m_e = q_w * p_e * q_d
m_d = q_w * q_e * p_d
m_we = p_w * p_e * q_d
m_wd = p_w * q_e * p_d
m_ed = q_w * p_e * p_d
m_wed = p_w * p_e * p_d

m_all_viburnum = c( m_0, m_w, m_e, m_d, m_we, m_wd, m_ed, m_wed )
names(m_all_viburnum) = c( "0", "w", "e", "d", "we", "wd", "ed", "wed")
print(m_all_viburnum)

# Network 2: rescaled viburnum
p_w = 8.95640194416046E-01
p_e = 3.22982251644135E-01
p_d = 9.7535240650177E-01

q_w = 1 - p_w
q_e = 1 - p_e
q_d = 1 - p_d

m_0 = q_w * q_e * q_d
m_w = p_w * q_e * q_d
m_e = q_w * p_e * q_d
m_d = q_w * q_e * p_d
m_we = p_w * p_e * q_d
m_wd = p_w * q_e * p_d
m_ed = q_w * p_e * p_d
m_wed = p_w * p_e * p_d

m_all_viburnum_2 = c( m_0, m_w, m_e, m_d, m_we, m_wd, m_ed, m_wed )
names(m_all_viburnum_2) = c( "0", "w", "e", "d", "we", "wd", "ed", "wed")
print(m_all_viburnum_2)
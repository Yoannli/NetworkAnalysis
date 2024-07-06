bg_color white
set grid_mode, 1
fetch 1UNR
fetch 2UZR
align (chain A & 1UNR), (chain A & 2UZR)
zoom (chain A & 1UNR)
remove resn HOH
set_color colordefault, [0.9657054978854287,0.9672433679354094,0.9680891964628989]
color colordefault, all
set_color RdBu_a, [0.02,0.188,0.38]
set_color RdBu_b, [0.127,0.396,0.669]
set_color RdBu_c, [0.263,0.576,0.765]
set_color RdBu_d, [0.566,0.769,0.869]
set_color RdBu_e, [0.82,0.898,0.941]
set_color RdBu_f, [0.966,0.967,0.968]
set_color RdBu_g, [0.992,0.859,0.78]
set_color RdBu_h, [0.955,0.642,0.506]
set_color RdBu_i, [0.839,0.376,0.302]
set_color RdBu_j, [0.692,0.092,0.168]
set_color RdBu_k, [0.404,0.0,0.122]
color colordefault, 1UNR & chain A
color colordefault, 2UZR & chain A
select a_grp_1_A, ((1UNR & chain A & resi 18,19,20,21,23,67,68,70))
color RdBu_a, a_grp_1_A
select b_grp_2_A, ((1UNR & chain A & resi 13,17,24,39,40,49,61,65,66,78,80,81,82,85,91,92,93,94,119,120) or (2UZR & chain A & resi 18,19,23,67,108,109,112,115,116,117,118,120))
color RdBu_b, b_grp_2_A
select c_grp_3_A, ((1UNR & chain A & resi 15,16,25,38,41,48,50,52,54,59,63,64,71,74,76,79,84,86,95,97,98,101,108,111,114,115,116,117,118) or (2UZR & chain A & resi 3,5,6,13,20,21,24,31,36,38,39,40,41,49,50,54,59,61,63,66,68,70,74,76,78,79,80,81,82,83,84,85,91,92,93,94,95,97,98,101,102,104,105,107,111,113,114,119))
color RdBu_c, c_grp_3_A
select d_grp_4_A, ((1UNR & chain A & resi 3,11,14,34,36,42,51,58,69,72,83,87,89,104,105,107,109,112,113) or (2UZR & chain A & resi 4,7,8,9,11,14,15,16,17,25,26,28,30,32,34,42,51,52,53,56,58,64,65,71,77,86))
color RdBu_d, d_grp_4_A
select e_grp_5_A, ((1UNR & chain A & resi 5,6,8,9,10,26,28,30,31,32,53,56,73,75,77,88,96,100) or (2UZR & chain A & resi 10,69,72,75,87,89,100,103,110))
color RdBu_e, e_grp_5_A
select f_grp_6_A, ((1UNR & chain A & resi 4,7,22,55,102,103,110) or (2UZR & chain A & resi 22,29,55,60,62,73,88,96,106))
color RdBu_f, f_grp_6_A
select g_grp_7_A, ((1UNR & chain A & resi 12,60,90,106) or (2UZR & chain A & resi 37,90,99))
color RdBu_g, g_grp_7_A
select h_grp_8_A, ((1UNR & chain A & resi 27,29,33,37,62,99) or (2UZR & chain A & resi 12,27))
color RdBu_h, h_grp_8_A
select i_grp_9_A, ((1UNR & chain A & resi 35) or (2UZR & chain A & resi 33,35))
color RdBu_i, i_grp_9_A
select k_grp_11_A, ((1UNR & chain A & resi 57) or (2UZR & chain A & resi 57))
color RdBu_k, k_grp_11_A
deselect

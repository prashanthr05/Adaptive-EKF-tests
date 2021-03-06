function h_withoutSkin = measurement_withoutSkin(v_Bx,v_By,v_Bz,omega_Bx,omega_By,omega_Bz,f_B_ox,f_B_oy,f_B_oz,mu_B_ox,mu_B_oy,mu_B_oz,f_B_cx,f_B_cy,f_B_cz,mu_B_cx,mu_B_cy,mu_B_cz,phi1,phi2,phi3,k_xx,k_yy,k_zz,c_xx,c_yy,c_zz,phi01,phi02,phi03,I_Bxx,I_Bxy,I_Bxz,I_Byy,I_Byz,I_Bzz,m,G_g1,G_g2,G_g3)
%MEASUREMENT_WITHOUTSKIN
%    H_WITHOUTSKIN = MEASUREMENT_WITHOUTSKIN(V_BX,V_BY,V_BZ,OMEGA_BX,OMEGA_BY,OMEGA_BZ,F_B_OX,F_B_OY,F_B_OZ,MU_B_OX,MU_B_OY,MU_B_OZ,F_B_CX,F_B_CY,F_B_CZ,MU_B_CX,MU_B_CY,MU_B_CZ,PHI1,PHI2,PHI3,K_XX,K_YY,K_ZZ,C_XX,C_YY,C_ZZ,PHI01,PHI02,PHI03,I_BXX,I_BXY,I_BXZ,I_BYY,I_BYZ,I_BZZ,M,G_G1,G_G2,G_G3)

%    This function was generated by the Symbolic Math Toolbox version 6.2.
%    23-Aug-2015 00:02:25

t2 = 1.0./m;
h_withoutSkin = [-f_B_cx.*t2+f_B_ox.*t2-omega_By.*v_Bz+omega_Bz.*v_By;-f_B_cy.*t2+f_B_oy.*t2+omega_Bx.*v_Bz-omega_Bz.*v_Bx;-f_B_cz.*t2+f_B_oz.*t2-omega_Bx.*v_By+omega_By.*v_Bx;omega_Bx;omega_By;omega_Bz;f_B_ox;f_B_oy;f_B_oz;mu_B_ox;mu_B_oy;mu_B_oz;f_B_cx;f_B_cy;f_B_cz;mu_B_cx;mu_B_cy;mu_B_cz];

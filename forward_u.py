# Define the weak form for forward PDE

def forward_pde_u(rho, u, v, f, s, bcs):
     # Solve for "u"
     a_forward_v = h_v(rho) * inner(sigma_v(u, Id), epsilon(v)) * dx
     a_forward_s = h_s(rho) * inner(sigma_s(u, Id), epsilon(v)) * dx
     a_forward_r = h_r(rho) * inner(sigma_r(u, Id), epsilon(v)) * dx
     a_forward = a_forward_v + a_forward_s + a_forward_r

     L_forward = inner(f, v) * ds(8) + s * h_r(rho) * inner(sigma_A(Id, Id), epsilon(v)) * dx
     L_forward_s = s * h_r(rho) * inner(sigma_A(Id, Id), epsilon(v)) * dx
     R_fwd = a_forward - L_forward
     R_fwd_s = a_forward - L_forward_s

     solve(R_fwd == 0, u, bcs = bcs)
     return u


def forward_pde_us(rho, u, v, f, s, bcs):
     # Solve for "u"
     a_forward_v = h_v(rho) * inner(sigma_v(u, Id), epsilon(v)) * dx
     a_forward_s = h_s(rho) * inner(sigma_s(u, Id), epsilon(v)) * dx
     a_forward_r = h_r(rho) * inner(sigma_r(u, Id), epsilon(v)) * dx
     a_forward = a_forward_v + a_forward_s + a_forward_r

     L_forward_s = s * h_r(rho) * inner(sigma_A(Id, Id), epsilon(v)) * dx
     R_fwd_s = a_forward - L_forward_s

     solve(R_fwd_s == 0, u, bcs = bcs)
     return u


def forward_pde_s(s, s_n, vs, dt, g, bcss):
     a_forward_s = s * vs * dx + dt * dot(grad(u), grad(vs)) * dx
     L_forward_s = (s_n + dt * g) * vs * dx

     R_heat_forward = a_forward_s - L_forward_s

     solve(R_heat_forward == 0, s, bcs = bcss)
     return s
def r_X(r_Y, e_XY, e_XY_prime):
    R_XY = (e_XY_prime / e_XY) - 1
    r_X_value = (r_Y - R_XY) / (1 + R_XY)
    return r_X_value
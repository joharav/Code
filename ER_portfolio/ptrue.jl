function ptrue(enn::Int64)
    pea = zeros(enn);
    pea[1] = 0.985         # beta
    pea[2] = 0.05          # delta
    pea[3] = 0.66          # rho_e
    pea[4] = 0.12          # sigma_e
    pea[5] =  0.587963      # nu
    pea[6] = 2.00          # gamma
    pea[7] = 0.059965      # f (durable adj cost share)
    pea[8] = 1             # w
    pea[9] = 0.0075        # r_foreign
    pea[10] = 5            # pd
    pea[11] = 0.722273    # kappa (asset adj cost)
    pea[12] = 0.0          # tau
    pea[13] = 1.0          # h
    pea[14] = 0.9          # rho_y
    pea[15] = 0.2          # sigma_y
    pea[16] = 0.531183           # chi (maintenance effectiveness)
    pea[17] = 0.397275         # ft  (time fixed cost)

    return pea
end

import numpy as np
import ot


def gwloss_partial(C1, C2, T, loss_fun='square_loss'):
    g = gwgrad_partial(C1, C2, T, loss_fun) * 0.5
    return np.sum(g * T)


def wloss(M, T):
    return np.sum(M * T)


def fgwloss_partial(alpha, M, C1, C2, T, loss_fun='square_loss'):
    return (1 - alpha) * wloss(M, T) + alpha * gwloss_partial(C1, C2, T, loss_fun)


def gwgrad_partial(C1, C2, T, loss_fun="square_loss"):
    """Compute the GW gradient, as one term in the FGW gradient.

    Note: we can not use the trick in Peyre16 as the marginals may not sum to 1.

    Parameters
    ----------
    C1: array of shape (n_p,n_p)
        intra-source cost matrix

    C2: array of shape (n_q,n_q)
        intra-target cost matrix

    T : array of shape(n_p, n_q)
        Transport matrix

    loss_fun

    Returns
    -------
    numpy.array of shape (n_p, n_q)
        gradient
    """
    if loss_fun == 'square_loss':
        def f1(a):
            return (a**2)

        def f2(b):
            return (b**2)

        def h1(a):
            return a

        def h2(b):
            return 2 * b
    elif loss_fun == 'kl_loss':
        def f1(a):
            return a * np.log(a + 1e-15) - a

        def f2(b):
            return b

        def h1(a):
            return a

        def h2(b):
            return np.log(b + 1e-15)

    #cC1 = np.dot(C1 ** 2 / 2, np.dot(T, np.ones(C2.shape[0]).reshape(-1, 1)))
    A = np.dot(
        f1(C1),
        np.dot(T, np.ones(C2.shape[0]).reshape(-1, 1))
    )

    #cC2 = np.dot(np.dot(np.ones(C1.shape[0]).reshape(1, -1), T), C2 ** 2 / 2)
    B = np.dot(
        np.dot(np.ones(C1.shape[0]).reshape(1, -1), T),
        f2(C2).T
    )  # TODO does f2(C2) here need transpose?

    constC = A + B
    #C = -np.dot(C1, T).dot(C2.T)
    C = -np.dot(h1(C1), T).dot(h2(C2).T)
    tens = constC + C
    return tens * 2


def fgwgrad_partial(alpha, M, C1, C2, T, loss_fun='square_loss'):
    return (1 - alpha) * M + alpha * gwgrad_partial(C1, C2, T, loss_fun)


def partial_fused_gromov_wasserstein(M, C1, C2, p, q, alpha, m=None, G0=None, loss_fun='square_loss', armijo=False, log=False, verbose=False, numItermax=1000, tol=1e-7, stopThr=1e-9, stopThr2=1e-9):
    if m is None:
        m = np.min((np.sum(p), np.sum(q)))
    elif m < 0:
        raise ValueError("Problem infeasible. Parameter m should be greater"
                         " than 0.")
    elif m > np.min((np.sum(p), np.sum(q))):
        raise ValueError("Problem infeasible. Parameter m should lower or"
                         " equal to min(|p|_1, |q|_1).")

    if G0 is None:
        G0 = np.outer(p, q)

    nb_dummies = 1
    dim_G_extended = (len(p) + nb_dummies, len(q) + nb_dummies)
    q_extended = np.append(q, [(np.sum(p) - m) / nb_dummies] * nb_dummies)
    p_extended = np.append(p, [(np.sum(q) - m) / nb_dummies] * nb_dummies)

    cpt = 0
    err = 1

    if log:
        log = {'err': [], 'loss': []}
    f_val = fgwloss_partial(alpha, M, C1, C2, G0, loss_fun)
    if verbose:
        print('{:5s}|{:12s}|{:8s}|{:8s}'.format(
            'It.', 'Loss', 'Relative loss', 'Absolute loss') + '\n' + '-' * 48)
        print('{:5d}|{:8e}|{:8e}|{:8e}'.format(cpt, f_val, 0, 0))

    while err > tol and cpt < numItermax:
        Gprev = np.copy(G0)
        old_fval = f_val

        gradF = fgwgrad_partial(alpha, M, C1, C2, G0, loss_fun)
        gradF_emd = np.zeros(dim_G_extended)
        gradF_emd[:len(p), :len(q)] = gradF
        gradF_emd[-nb_dummies:, -nb_dummies:] = np.max(gradF) * 1e2
        gradF_emd = np.asarray(gradF_emd, dtype=np.float64)

        Gc, logemd = ot.lp.emd(p_extended, q_extended, gradF_emd, log=True)
        if logemd['warning'] is not None:
            raise ValueError("Error in the EMD resolution: try to increase the"
                             " number of dummy points")

        G0 = Gc[:len(p), :len(q)]

        if cpt % 10 == 0:  # to speed up the computations
            err = np.linalg.norm(G0 - Gprev)
            if log:
                log['err'].append(err)
            # if verbose:
            #     if cpt % 200 == 0:
            #         print('{:5s}|{:12s}|{:12s}'.format(
            #             'It.', 'Err', 'Loss') + '\n' + '-' * 31)
            #         print('{:5d}|{:8e}|{:8e}'.format(cpt, err,
            #                                          fgwloss_partial(alpha, M, C1, C2, G0, loss_fun)))

        deltaG = G0 - Gprev

        if not armijo:
            a = alpha * gwloss_partial(C1, C2, deltaG, loss_fun)
            b = (1 - alpha) * wloss(M, deltaG) + 2 * alpha * np.sum(gwgrad_partial(C1, C2, deltaG, loss_fun) * 0.5, Gprev)
            # c = (1 - alpha) * wloss(M, Gprev) + alpha * gwloss_partial(C1, C2, Gprev, loss_fun)
            c = fgwloss_partial(alpha, M, C1, C2, Gprev, loss_fun)

            gamma = ot.optim.solve_1d_linesearch_quad(a, b, c)
            # f_val = a * gamma ** 2 + b * gamma + c
        else:
            def f(x, alpha, M, C1, C2, lossfunc):
                return fgwloss_partial(alpha, M, C1, C2, x, lossfunc)
            xk = Gprev
            pk = deltaG
            gfk = fgwgrad_partial(alpha, M, C1, C2, xk, loss_fun)
            old_val = fgwloss_partial(alpha, M, C1, C2, xk, loss_fun)
            args = (alpha, M, C1, C2, loss_fun)
            gamma, fc, fa = ot.optim.line_search_armijo(f, xk, pk, gfk, old_val, args)
            # f_val = f(xk + gamma * pk, alpha, M, C1, C2, loss_fun)

        if gamma == 0:
            cpt = numItermax
        G0 = Gprev + gamma * deltaG
        f_val = fgwloss_partial(alpha, M, C1, C2, G0, loss_fun)
        cpt += 1

        # TODO: better stopping criteria?
        abs_delta_fval = abs(f_val - old_fval)
        relative_delta_fval = abs_delta_fval / abs(f_val)
        if relative_delta_fval < stopThr or abs_delta_fval < stopThr2:
            cpt = numItermax
        if log:
            log['loss'].append(f_val)
        if verbose:
            if cpt % 20 == 0:
                print('{:5s}|{:12s}|{:8s}|{:8s}'.format(
                    'It.', 'Loss', 'Relative loss', 'Absolute loss') + '\n' + '-' * 48)
            print('{:5d}|{:8e}|{:8e}|{:8e}'.format(cpt, f_val, relative_delta_fval, abs_delta_fval))

    if log:
        log['partial_fgw_cost'] = fgwloss_partial(alpha, M, C1, C2, G0, loss_fun)
        return G0[:len(p), :len(q)], log
    else:
        return G0[:len(p), :len(q)]










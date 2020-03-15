import numpy as np

BACKTRACK_EPS = 1E-3


def backtrack(target_function, fval, old_fval, new_gradient,
              step_direction, armijo_line_search_factor,
              step_size, step_size_adjustment, reset_step_size,
              debugging, ALPHA_CONST, RESET_CONST, accelerate,
              linesearch):
    '''
    A wrapper for the backtracking algorithms employed in CG, SD,
    BFGS, and LBFGS.  Not intended for use beyond these codes.
    '''

    # Get function to describe linesearch
    if linesearch is 'armijo':
        def check_backtrack(f1,
                            f0,
                            gk,
                            pk,
                            armijio_line_search_factor,
                            step_size):
            a = f1 - f0
            b = armijio_line_search_factor * step_size * np.dot(gk, pk)
            return a > b
    elif linesearch is 'backtrack':
        def check_backtrack(f1,
                            f0,
                            pk,
                            gk,
                            armijio_line_search_factor,
                            step_size):
            a = f1 - f0
            a /= (abs(f1) + abs(f0))
            b = BACKTRACK_EPS
            return a > b
    else:
        def check_backtrack(f1,
                            f0,
                            pk,
                            gk,
                            armijio_line_search_factor,
                            step_size):
            return False

    cont_flag = False
    # Here we deal with backtracking if needed
    if (target_function is not None and
        check_backtrack(fval,
                        old_fval,
                        new_gradient,
                        step_direction,
                        armijo_line_search_factor,
                        step_size)):

        step_size *= np.float64(step_size_adjustment)

        # Reset the Inverse Hessian if desired.
        # It is still up for debate if this is to be recommended or not.
        # As the inverse hessian corects itself, it might not be important
        # to do this.
        reset_step_size = RESET_CONST
        if debugging:
            print("\tBACKTRACKING! Step Size is %lg" % step_size)
        cont_flag = True

    # This allows for the edge case in which after decreasing step_size,
    # a situation arises in which larger alphas are acceptable again.
    # Thus, we reset to the original step_size
    elif reset_step_size is not None:
        reset_step_size -= 1
        # If we want to reset_step_size and step_size has been decreased
        # before, set to initial vals
        if reset_step_size < 0 and step_size < ALPHA_CONST:
            step_size, reset_step_size = ALPHA_CONST, RESET_CONST
            # Once again, debatable if we want this here.  When reseting
            # step sizes we might have a better H inverse than the Identity
            # would be.
            if debugging:
                    print("\tACCELERATING 1! Step Size is %lg" % step_size)
            cont_flag = True

        # If we want to reset_step_size and we've never decreased before,
        # we can take larger steps. We increase step sizes by a factor
        # of 1/step_size_adjustment
        elif (reset_step_size < 0 and
              step_size >= ALPHA_CONST and
              accelerate):
            if linesearch is not None:
                if debugging:
                    print("\tACCELERATING 2! Step Size is %lg" % step_size)
                step_size /= step_size_adjustment
            reset_step_size = RESET_CONST

    return reset_step_size, step_size, cont_flag

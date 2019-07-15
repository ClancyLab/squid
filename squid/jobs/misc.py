def close_pipes(p):
    """
    A simple function to close the pipes if they remain open.
    """
    if p is not None:
        if p.stdout is not None:
            p.stdout.close()
        if p.stderr is not None:
            p.stderr.close()

import time


class JobObject(object):
    '''
    Job class to wrap simulations for queue submission.

    **Parameters**

        name: *str*
            Name of the simulation on the queue.
        process_handle: *process_handle, optional*
            The process handle, returned by subprocess.Popen.
        job_id: *str, optional*
            The job id.  Usually this should be unique.

    **Returns**

        job_obj: :class:`squid.jobs.container.JobObject`
            A Job object.
    '''

    def __init__(self, name, process_handle=None, job_id=None):
        self.name = name
        self.process_handle = process_handle
        self.job_id = job_id

    def wait(self, tsleep=60, verbose=False):
        '''
        Hang until simulation has finished.

        tsleep: *int, optional*
            How long to wait before checking if the job has finished in
            the loop.  Default is 1 minute.
        verbose: *bool, optional*
            Whether to print repeatedly on each check or not.

        **Returns**

            None
        '''
        if self.process_handle is not None:
            self.process_handle.wait()
        else:
            while True:
                if not self.is_finished():
                    if verbose:
                        print("Job (%s) is still running..." % self.name)
                    time.sleep(tsleep)
                else:
                    break

    def get_all_jobs(self, detail=3):
        '''
        Get a list of all jobs that are running and/or pending.

        **Parameters**

            detail: *int, optional*
                How much detail to get when finding jobs on the queue.

        **Returns**

            jobs_on_queue: *list, ...*
                A list of all jobs on the queue, and any other relevant
                information requested.
        '''
        return []

    def is_finished(self):
        '''
        Check if simulation has finished or not.

        **Returns**

            is_on_queue: *bool*
                Whether the simulation is still running (True), or not (False).
        '''
        if self.process_handle is not None:
            return self.process_handle.poll() == 0
        if self.job_id is not None:
            running = any([
                self.job_id in j for j in self.get_all_jobs(detail=3)])
        else:
            running = self.name in self.get_all_jobs(detail=0)
        return not running

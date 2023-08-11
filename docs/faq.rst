.. _faq::

Frequently-Asked Questions
==========================

I got this error/warning message...
-----------------------------------


   **WARNING**: LOCAL JWST PRD VERSION PRDOPSSOC-059 CANNOT BE CHECKED AGAINST ONLINE VERSION


Likely not connected to the internet or change in `pysiaf <https://github.com/spacetelescope/pysiaf>`_.  Verify the connectivity and/or upgrade the ``pysiaf`` version.


    RuntimeError: 
    An attempt has been made to start a new process before the
    current process has finished its bootstrapping phase.

    This probably means that you are not using fork to start your
    child processes and you have forgotten to use the proper idiom
    in the main module:
            if __name__ == '__main__':
                freeze_support()


    The "freeze_support()" line can be omitted if the program
    is not going to be frozen to produce an executable.


Likely running in the __main__ of a script. Try putting commands inside of a code block, such as:

.. code:: python
    if __name__=='__main__':  # doctest: +SKIP
        my_commands_here()  # doctest: +SKIP


Why do I....
------------

* **Why do I get the same logging message printed multiple times?**  This is likely because there are multiple instances of the logger running for the separate threads.
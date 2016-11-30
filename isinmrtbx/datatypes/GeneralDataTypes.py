from __future__ import print_function

AUTHOR = 'Tomas Psorn'
CONTACT = 'tomaspsorn@isibrno.cz'


class ParameterClass():

    def __init__(self):
        pass

    def export(self, show=True):
        """

        Returns
        -------
        if show:
            prints all parameters and their values
        else:
            dictionary of parameters
        """

        dict = {}

        for key in dir(self):
            dict[key] = eval("self." + key)

            if show:
                print(key)
                print(dict[key])
                print()

        if show:
            return
        else:
            return dict


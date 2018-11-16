from __future__ import absolute_import, print_function, division

import numpy as np
import astropy.units as u

from pint.models.timing_model import DelayComponent, MissingParameter
from pint.models import parameter as p


class IFunc(DelayComponent):
    """This class provides timing noise using tabulated IFUNC values."""
    register = True
    category = 'ifunc'
    def __init__(self):
        super(IFunc, self).__init__()

        self.add_param(p.prefixParameter(name="IFUNC1", units='s',
                                         description="IFUNC components",
                                         type_match='pair', long_double=True,
                                         parameter_type='pair'))
        self.delay_funcs_component += [self.ifunc_delay,]
        self.category = 'ifunc'

    def setup(self):
        super(IFunc, self).setup()
        ifunc_terms = list(self.get_prefix_mapping_component('IFUNC').keys())
        ifunc_terms.sort()
        ifunc_in_order = list(range(1, max(ifunc_terms)+1))
        if not ifunc_terms == ifunc_in_order:
            diff = list(set(ifunc_in_order) - set(ifunc_terms))
            raise MissingParameter("IFunc", "IFUNC%d"%diff[0])

        self.num_ifunc_terms = len(ifunc_terms)

    def print_par(self, ):
        result = ''
        ifunc_terms = ["IFUNC%d" % ii for ii in
                      range(1, self.num_ifunc_terms + 1)]

        for ft in ifunc_terms:
            par = getattr(self, ft)
            result += par.as_parfile_line()

        return result

    def ifunc_delay(self, toas, acc_delay=None):
        delays = 0
        ifunc_names = ["IFUNC%d" % ii for ii in
                      range(1, self.num_ifunc_terms + 1)]
        ifunc_terms = [getattr(self, name) for name in ifunc_names]
        time = self.barycentric_time = toas.get_mjds().value

        ifunc_mjd = []
        ifunc_v = []
        for ft in ifunc_terms:
            mjd, v = ft.quantity
            ifunc_mjd.append(float(mjd.value))
            ifunc_v.append(float(v.value))
        ifunc_mjd = np.array(ifunc_mjd)
        ifunc_v = np.array(ifunc_v)

        delays = -np.interp(time, ifunc_mjd, ifunc_v) * u.s
        return delays

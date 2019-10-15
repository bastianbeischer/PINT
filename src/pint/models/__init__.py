"""Implementations of pulsar timing models."""
# Import the main timing model classes
from pint.models.timing_model import TimingModel

# Import all standard model components here
from pint.models.astrometry import AstrometryEquatorial, AstrometryEcliptic
from pint.models.binary_bt import BinaryBT
from pint.models.binary_dd import BinaryDD
from pint.models.binary_ell1 import BinaryELL1, BinaryELL1H
from pint.models.binary_ddk import BinaryDDK
from pint.models.dispersion_model import DispersionDM, DispersionDMX
from pint.models.solar_wind_dispersion import SolarWindDispersion
from pint.models.spindown import Spindown
from pint.models.frequency_dependent import FD
from pint.models.absolute_phase import AbsPhase
from pint.models.glitch import Glitch
from pint.models.jump import DelayJump, PhaseJump
from pint.models.solar_system_shapiro import SolarSystemShapiro
from pint.models.noise_model import ScaleToaError, EcorrNoise, PLRedNoise
from pint.models.model_builder import get_model
from pint.models.wave import Wave
from pint.models.ifunc import IFunc

# Define a standard basic model
StandardTimingModel = TimingModel("StandardTimingModel",
          (AstrometryEquatorial(), Spindown(), DispersionDM(), SolarSystemShapiro()))
# BTTimingModel = generate_timing_model("BTTimingModel",
#         (Astrometry, Spindown, Dispersion, SolarSystemShapiro, BT))
# DDTimingModel = generate_timing_model("DDTimingModel",
#         (Astrometry, Spindown, Dispersion, SolarSystemShapiro, DD))

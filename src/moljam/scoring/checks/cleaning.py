from .cleaning_core import CleaningCoreMixin
from .cleaning_save import CleaningSaveMixin


class CleaningMixin(
    CleaningCoreMixin,
    CleaningSaveMixin,
):
    pass


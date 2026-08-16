from .._common import *


class ExperimentalMethodChecksMixin:
    def check_experimental_methods(self):
        """Compatibility no-op: experimental method analysis has been removed."""
        print("Experimental method analysis is disabled and no longer recorded")
        self.completed_checks.add('check_experimental_methods')
        return 0

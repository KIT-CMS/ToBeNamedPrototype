import logging

log = logging.getLogger(__name__)


class Quantity:
    def __init__(self, name):
        self.name = name
        # structure for storing shifts:
        # {scope1: {shift1 : name2,
        #          shift2: name2},
        #  scope2: {shift1 : name2,
        #           shift2: name2}, ...}
        self.shifts = {}
        self.ignored_shifts = {}
        self.children = {}
        self.defined_for_scopes = []
        log.debug("Setting up new Quantity {}".format(self.name))

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name

    def reserve_scope(self, scope):
        """
        Function to reserve a scope for a given quantity. The scopes, in which a quantity is used
        as an output are tracked in the output_scopes list.
        If a quantity is already used within a given scope as output, this will result in an exception.
        This check is triggered for every Producer.
        """
        log.debug("Checking {} / scope {}".format(self.name, scope))
        if scope not in self.defined_for_scopes:
            self.defined_for_scopes.append(scope)
        else:
            log.error(
                "Quantity {} is already defined in {} scope !".format(self.name, scope)
            )
            raise Exception

    def get_leaf(self, shift, scope):
        """
        Function to get the leaf of a given shift within a given scope.
        A leaf is the name of the quantity used for that scope/shift combination.
        If no shift is defined, the name of the quantity is returned.

        Args:
            shift (str): Name of the shift for which the leaf should be returned
            scope (str): Scope for which the leaf should be returned
        Returns:
            str. Name of the leaf
        """
        log.debug("{} - getting shift {} for scope {}".format(self.name, shift, scope))
        leaf = self.name
        if scope in self.shifts.keys():
            if shift in self.shifts[scope].keys():
                leaf = self.shifts[scope][shift]
        return leaf

    def get_leafs_of_scope(self, scope):
        """
        Function returns a list of all leafs, which are defined for a given scope.

        Args:
            scope (str): Scope for which the leafs should be returned
        Returns:
            list. List of leafs
        """
        result = [self.name] + [
            self.get_leaf(shift, scope) for shift in self.get_shifts(scope)
        ]
        return result

    def shift(self, name, scope):
        """
        Function to define a shift for a given scope. If the shift is marked as ignored, nothing will be added.
        When a new shift is defined, all child quantities of a given quantity will be shifted as well.
        Shifts defined in the global scope, are added to the shift dictionary of all scopes.
        Shifts that are exclusive to a single scope are only added to the given scope.

        Args:
            name (str): Name of the shift
            scope (str): Scope for which the shift should be defined
        Returns:
            None
        """
        if scope in self.ignored_shifts.keys():
            if name in self.ignored_shifts[scope]:
                log.debug("Ignoring shift {} for quantity {}".format(name, self.name))
                return
        log.debug("Adding shift {} to quantity {}".format(name, self.name))
        # adding new shifts to scopes:
        # if a shift is defined for the global scope, it should also be added for all other scopes,
        # therefore, if a new scope is added,
        # this scope is a copy of the global scope.
        if scope not in self.shifts.keys():
            if "global" in self.shifts.keys():
                # make a copy of the global shifts
                self.shifts[scope] = self.shifts["global"]
            else:
                self.shifts[scope] = {}
        scopes = [scope]
        if scope == "global" and scope in self.shifts.keys():
            # in this case, we add the shift to all existing scopes
            scopes = self.shifts.keys()
        for scope in scopes:
            if name not in self.shifts[scope]:
                self.shifts[scope][name] = self.name + name
                if scope == "global":  # shift children in all scopes if scope is global
                    for any_scope in self.children:
                        for c in self.children[any_scope]:
                            c.shift(name, any_scope)
                else:
                    if scope in self.children.keys():
                        for c in self.children[scope]:
                            c.shift(name, scope)

    def ignore_shift(self, name, scope):
        """
        Function to ignore a shift for a given scope.

        Args:
            name (str): Name of the shift to be ignored
            scope (str): Scope for which the shift should be ignored
        Returns:
            None
        """
        log.debug("Make quantity {} ignore shift {}".format(self.name, name))
        if scope not in self.ignored_shifts.keys():
            self.ignored_shifts[scope] = {}
        self.ignored_shifts[scope] = name

    def copy(self, name):
        """
        Generate a copy of the current quantity with a new name.

        Args:
            name (str): Name of the new quantity
        Returns:
            Quantity. a new Quantity object.
        """
        copy = Quantity(name)
        copy.shifts = self.shifts
        copy.children = self.children
        copy.ignored_shifts = self.ignored_shifts
        copy.defined_for_scopes = self.defined_for_scopes
        return copy

    def adopt(self, child, scope):
        """
        Adopt a child quantity to the current quantity.
        An adopted quantity will inherit all shifts of the partent quantity.

        Args:
            child (Quantity): The child quantity
            scope (str): Scope for which the child should be adopted
        Returns:
            None
        """
        log.debug(
            "Adopting child quantity {} to quantity {} in scope {}".format(
                child.name, self.name, scope
            )
        )
        if scope not in self.children.keys():
            self.children[scope] = []
        self.children[scope].append(child)

    def get_shifts(self, scope):
        """
        Function returns a list of all shifts, which are defined for a given scope.

        Args:
            scope (str): Scope for which shifts should be returned
        Returns:
            list: List of all shifts, which are defined for a given scope.
        """
        if scope in self.shifts.keys():
            return list(self.shifts[scope].keys())
        else:
            return []


class QuantityGroup(Quantity):
    # A Quantity Group is a group of quantities, that all have the same settings, but different names.
    def __init__(self, name):
        super().__init__(name)
        self.quantities = []

    def copy(self, name):
        log.error("Copy is not allowed for a Quantity Group !")
        raise Exception

    def add(self, name):
        # add a new Quantity to the group. This quantity contains the identical shifts as the group itself
        quantity = Quantity(name)
        quantity.shifts = self.shifts
        quantity.children = self.children
        quantity.ignored_shifts = self.ignored_shifts
        quantity.defined_for_scopes = self.defined_for_scopes
        self.quantities.append(quantity)

    def get_leafs_of_scope(self, scope):
        # For the writeout, we have to loop over all quantities in the group and return them all (plus their shifts) in a list
        output = []
        for quantity in self.quantities:
            output.extend(quantity.get_leafs_of_scope(scope))
        return output


class NanoAODQuantity(Quantity):
    def __init__(self, name):
        super().__init__(name)

    # Quantities from the NanoAOD are not designed to be directly usable as output
    def reserve_scope(self, scope):
        log.error(
            "Quantity {} is a NanoAOD quantity and cant be used as output !".format(
                self.name
            )
        )
        raise Exception

    # if the shifted version of a quantity already exists in the input,
    # this function can be used to register an
    # branch from the input as a shifted version of a quantity
    def register_external_shift(self, shift, scope, external_name):
        if scope not in self.shifts.keys():
            self.shifts[scope] = {}
        if scope == "global":
            for allscope in self.shifts.keys():
                self.shifts[allscope][shift] = external_name
        else:
            self.shifts[scope][shift] = external_name

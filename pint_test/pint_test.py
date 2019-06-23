from pint import UnitRegistry

#%%
ureg = UnitRegistry()

distance = 24.0 * ureg.meter
print(distance)

time = 8.0 * ureg.second

print(repr(time))

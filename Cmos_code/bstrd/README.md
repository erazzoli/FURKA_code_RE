# BStrd

Online analysis via bsread

## Examples

### Creating channels

```python
from bstrd import BS, bsstream

pid = BS("pid")
trace = BS("SATES21-GES1:A1_VALUES")
```

with optional `modulo` and `offset`

```python
inten = BS("SATFE10-PEPG046:FCUP-INTENSITY-CAL", modulo=10, offset=1)
```

### Receiving data

Read from channel:

```python
for _ in bsstream:
    val = trace.value
    print(val)
```


Read from data dict:

```python
for data in bsstream:
    val = data["SATES21-GES1:A1_VALUES"]
    print(val)
```


Iterate a single channel:

```python
for val in trace:
    print(val)
```



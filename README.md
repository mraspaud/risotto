# risotto
Rice aka USES aka SZIP decoder in Python/Cython

This an implementation of CCSDS's document 121.0-B-2 "Lossless Data Compression" (Blue Book)
 http://public.ccsds.org/publications/archive/121x0b2.pdf

Usage:

```python
from risotto import decode

decoded_data = decode(data, N, J, S)
```

where:
  N is the bits per sample size (eg 15)
  J is the block size (eg 8)
  S ih the segment size (eg 128)




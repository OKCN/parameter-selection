Parameter Scripts
===========================

Files
-----

Most of files are modified or copied from project [lwe-frodo](https://github.com/lwe-frodo/parameter-selection).

The main standalone scripts in this directory are:

OKCN-LWR
- **OKCN-LWR.cpp** - error rate of OKCN-LWR (Table 2)
- **LWR-SEC.py** - concrete security estimation of LWR (Table 3)

OKCN-LWE
- **OKCN-LWE.py** - OKCN-LWE (Table 6 and Table 7)
- **OKCN-LWE-512.py** - OKCN-LWE 512-bit (last line of Table 6)
- **Cut-LWE.py** - OKCN-LWE with cutting some least significant bits (Table 8)

OKCN/AKCN-RLWE
- **OKCN-AKCN-RLWE-CODE.py** - per bit error rate of OKCN/AKCN-RLWE, and variant with SEC code (Table 9)

Compile
-------

Use following command to compile the OKCN-LWR.cpp
```Bash
g++ -O2 OKCN-LWR.cpp -o OKCN-LWR -Wall
```
It may take you some minutes to run this program.


License
-------

This software is licensed under the MIT License.  
For details, see [LICENSE.txt](https://github.com/lwe-frodo/parameter-selection/blob/master/LICENSE.txt).

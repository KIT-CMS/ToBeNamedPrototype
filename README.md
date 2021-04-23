# ToBeNamedPrototype

## Environment setup

Current requirement is a working mount to CERN's `/cvmfs/sft.cern.ch`

```bash
source init.sh
```

## Compilation for main executable

```bash
g++ analysis.cxx $(root-config --cflags --libs) -fconcepts
```

## Profiling with perf & flamegraph for CPU

Running profiling on executable

```bash
perf record ./a.out
```

Print out the report

```bash
perf report
```

Get flamegraph repo

```bash
git clone https://github.com/brendangregg/FlameGraph
```

Create flamegraph

```bash
perf script > out.perf
FlameGraph/stackcollapse-perf.pl out.perf > out.folded
FlameGraph/flamegraph.pl out.folded > out.svg
```

## Profiling with valgrind massif for Memory

```bash
valgrind --tool=massif ./a.out
ms_print massif.out.4103388 > massif.log
```

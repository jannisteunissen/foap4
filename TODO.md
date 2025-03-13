# Ideas to implement

## Simplifications

- Do not split quadrant families between ranks

## Indexing quadrants

On the GPU, use a fixed linear array. p4est dynamically changes quadrant array when coarsening/refining or load balancing so that it does not contain holes.

Simple strategy: use same indices on GPU as in p4est (global index, cumulative over trees). Copy all data after a change, which is inefficient. Use two copies of GPU array, one with multiple time states, one with a single state.

## Mesh adaptation

Employ same approach as in userdata examples

## Partitioning

Can use something like p4est_transfer_fixed, but (in future) replace MPI calls with GPU ones.

## Ghost cell exchange

Use p4est_iterate with the face callback.

1. Count the number of faces of different kinds (same level, coarse-to-fine, fine-to-coarse).
2. For each type of face, store how many to send/recv to each other proc

Then, run with a different callback function that populates arrays with communication info and local exchanges.

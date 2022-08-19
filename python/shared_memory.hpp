#include <stdio.h>
#include <stdlib.h>
#include <iostream>


#ifdef HAVE_MPI

  #include <mpi.h>
  void MPI_configuration(MPI_Comm *allcomm, MPI_Comm *nodecomm, int*node_size, int*node_rank, int*world_size, int*world_rank)
      {
        // Read MPI properties
        *allcomm = MPI_COMM_WORLD;
        int MPI_initialized_already;

        MPI_Initialized(&MPI_initialized_already);
        if (!MPI_initialized_already)
            MPI_Init(NULL, NULL);
        

        MPI_Comm_size(*allcomm, world_size);
        MPI_Comm_rank(*allcomm, world_rank);

        // Create node-local communicator
        MPI_Comm_split_type(*allcomm, MPI_COMM_TYPE_SHARED, *world_rank, MPI_INFO_NULL, nodecomm);
        MPI_Comm_size(*nodecomm, node_size);
        MPI_Comm_rank(*nodecomm, node_rank);
    }

  template<class T>
  T* SharedMemoryArray(size_t arraysize) {
      // Array properties
      T *array, *localarray;
      int localarraysize;
      
      // MPI configuration:
      MPI_Comm allcomm, nodecomm;
      int node_size, node_rank;
      int world_size, world_rank;
      MPI_configuration(&allcomm, &nodecomm, &node_size, &node_rank, &world_size, &world_rank);

      // Shared memmory vars
      int *model; 
      int i, flag;
      int windisp;
      int *winptr;
      MPI_Aint winsize;
      MPI_Win winarray;
      
      // Only rank 0 on a node actually allocates memory
      localarraysize = 0;
      if (node_rank == 0) localarraysize = arraysize;

      MPI_Win_allocate_shared(localarraysize*sizeof(int), sizeof(int),
              MPI_INFO_NULL, nodecomm, &localarray, &winarray);

      MPI_Win_get_attr(winarray, MPI_WIN_MODEL, &model, &flag);

      // need to get local pointer valid for array on rank 0

      array = localarray;

      if (node_rank != 0)
      {
          MPI_Win_shared_query(winarray, 0, &winsize, &windisp, &array);
      }

      MPI_Win_fence(0, winarray);

      return array;
          
    }
  
  template<class T>
  void memcpy_shared(T* dest, const T* src, size_t count ) 
  {
    // MPI configuration:
    MPI_Comm allcomm, nodecomm;
    int node_size, node_rank;
    int world_size, world_rank;
    MPI_configuration(&allcomm, &nodecomm, &node_size, &node_rank, &world_size, &world_rank);

    if (node_rank == 0) {
        memcpy(dest, src, count);
        std::cout << "copying the values from process ["<<node_rank << "/" <<node_size<<"] \n" ;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
#else
    
  template<class T>
  T* SharedMemoryArray(size_t arraysize) {
    return new T[arraysize];
    }
  
  template<class T>
  void memcpy_shared(T* dest, const T* src, size_t count ) {
    std::cout << "copying the values from all processes \n" ;
    memcpy(dest, src, count);
    }

#endif











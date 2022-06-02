#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkCPPythonAdaptorAPI.h"

#include "vtkDoubleArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkSmartPointer.h"
#include "vtkUniformGrid.h"
#include "vtkRectilinearGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPointData.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkMultiPieceDataSet.h"
#include "vtkCellType.h"
#include "vtkCellData.h"
#include "vtkSOADataArrayTemplate.h"
#include "vtkNew.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <array>

#ifdef __cplusplus
extern "C" {
#endif
  void add_input_description(const char* grid_name);
  void set_time_data(int* step, double* time);
  void build_uniform_grid(int* rank, int* size, int* extent, int* zero_extent, int* whole_extent, double* origin, double* spacing, const char* grid_name);
  void build_rectilinear_grid(int* rank, int* size, int* extent, double* x_coords, double* y_coords, double* z_coords, int* ghost_levels, const char* grid_name);
  void build_unstructured_grid(long long int* count, const char* grid_name);
  void set_scalar_field(int* rank, double* scalar, char* name, const char* grid_name);
  void set_vector_field(int* rank, double* x_comp, double* y_comp, double* z_comp, char* name, const char* grid_name);
  void set_particle_position(long long int* id, double* x, double* y, double* z, const char* grid_name);
  void set_scalar_particle(long long int* id, double* scalar, char* name, const char* grid_name);
  void set_vector_particle(long long int* id, double* x_comp, double* y_comp, double* z_comp, char* name, const char* grid_name);
#ifdef __cplusplus
}
#endif

void add_input_description(const char* grid_name) {   
  if(!vtkCPPythonAdaptorAPI::GetCoProcessorData()) {
    vtkGenericWarningMacro("unable to access coprocessor data");
    return;
  }
  if(!vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName(grid_name)) {
    vtkCPPythonAdaptorAPI::GetCoProcessorData()->AddInput(grid_name);
  }
}

void set_time_data(int* step, double* time) {
  if(!vtkCPPythonAdaptorAPI::GetCoProcessorData()) {
    vtkGenericWarningMacro("unable to access coprocessor data");
    return;
  }
  vtkCPPythonAdaptorAPI::GetCoProcessorData()->SetTimeData(*time, static_cast<vtkIdType>(*step));
}

void build_uniform_grid(int* rank, int* size, int* extent, int* zero_extent, int* whole_extent, double* origin, double* spacing, const char* grid_name) {
  if(!vtkCPPythonAdaptorAPI::GetCoProcessorData()) {
    vtkGenericWarningMacro("unable to access coprocessor data");
    return;
  }
  vtkSmartPointer<vtkCPInputDataDescription> idd = vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName(grid_name);
  if(!idd) {
    vtkGenericWarningMacro("cannot access data description to attach grid to");
    return;
  }
  if(!idd->GetGrid()) {
    vtkSmartPointer<vtkImageData> local_grid = vtkSmartPointer<vtkImageData>::New();
    local_grid->SetOrigin(origin);
    local_grid->SetSpacing(spacing);
    local_grid->SetExtent(extent);
    //local_grid->GenerateGhostArray(zero_extent);
    idd->SetWholeExtent(whole_extent);
    idd->SetGrid(local_grid);
  } 
}

void build_rectilinear_grid(int* rank, int* size, int* extent, double* x_coords, double* y_coords, double* z_coords, int* ghost_levels, const char* grid_name) {  
  if(!vtkCPPythonAdaptorAPI::GetCoProcessorData()) {
    vtkGenericWarningMacro("unable to access coprocessor data");
    return;
  }
  vtkSmartPointer<vtkCPInputDataDescription> idd = vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName(grid_name);
  if(!idd) {
    vtkGenericWarningMacro("cannot access data description to attach grid to");
    return;
  }
  if(!idd->GetGrid()) {
    vtkSmartPointer<vtkRectilinearGrid> local_grid = vtkSmartPointer<vtkRectilinearGrid>::New();
    local_grid->SetExtent(extent);
    vtkDoubleArray* x_array = vtkDoubleArray::New();
    vtkDoubleArray* y_array = vtkDoubleArray::New();
    vtkDoubleArray* z_array = vtkDoubleArray::New();
    x_array->SetNumberOfComponents(1);
    y_array->SetNumberOfComponents(1);
    z_array->SetNumberOfComponents(1);
    x_array->SetVoidArray(x_coords, static_cast<vtkIdType>(extent[1] - extent[0] + 1), 1);
    y_array->SetVoidArray(y_coords, static_cast<vtkIdType>(extent[3] - extent[2] + 1), 1);
    z_array->SetVoidArray(z_coords, static_cast<vtkIdType>(extent[5] - extent[4] + 1), 1);
    local_grid->SetXCoordinates(x_array);
    local_grid->SetYCoordinates(y_array);
    local_grid->SetZCoordinates(z_array);
    x_array->Delete();
    y_array->Delete();
    z_array->Delete();
    
    //vtkUnsignedCharArray* ghost_cells = vtkUnsignedCharArray::New();
    //ghost_cells->SetNumberOfTuples(local_grid->GetNumberOfCells());
    //ghost_cells->SetName(vtkDataSetAttributes::GhostArrayName());
    //ghost_cells->Fill(0);
    //vtkIdType counter = 0;
    //for(int k = 0; k < extent[5] - extent[4]; k++) {
    //  for(int j = 0; j < extent[3] - extent[2]; j++) {
    //    for(int i = 0; i < extent[1] - extent[0]; i++) {
    //      if(k < *ghost_levels || k > (extent[5] - extent[4] - *ghost_levels) || j < *ghost_levels || j > (extent[3] - extent[2] - *ghost_levels) || i < *ghost_levels || i > (extent[1] - extent[0] - *ghost_levels)) {
    //        unsigned char value = vtkDataSetAttributes::DUPLICATECELL;
    //        ghost_cells->SetTypedTuple(counter, &value);
    //      }
    //      counter++;
    //    }
    //  }
    //}
    //local_grid->GetCellData()->AddArray(ghost_cells);
    //ghost_cells->Delete();
    
    vtkSmartPointer<vtkMultiPieceDataSet> multi_piece = vtkSmartPointer<vtkMultiPieceDataSet>::New();
    multi_piece->SetNumberOfPieces(*size);
    multi_piece->SetPiece(*rank, local_grid);
    //local_grid->Delete();
    vtkSmartPointer<vtkMultiBlockDataSet> grid = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    grid->SetNumberOfBlocks(1);
    grid->SetBlock(0, multi_piece);
    //multi_piece->Delete();
    idd->SetGrid(grid);
    //grid->Delete();
  } else {
  //
  }
}

void build_unstructured_grid(long long int* count, const char* grid_name) {
  if(!vtkCPPythonAdaptorAPI::GetCoProcessorData()) {
    vtkGenericWarningMacro("unable to access coprocessor data");
    return;
  }
  vtkSmartPointer<vtkCPInputDataDescription> idd = vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName(grid_name);
  if(!idd) {
    vtkGenericWarningMacro("cannot access data description");
    return;
  }
  if(!idd->GetGrid()) {
    vtkPoints* points = vtkPoints::New();
    points->SetNumberOfPoints(static_cast<vtkIdType>(*count));
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);
    points->Delete();
    grid->Allocate(grid->GetNumberOfPoints());
    vtkIdType cell_id[1];
    for(vtkIdType i = 0; i < static_cast<vtkIdType>(*count); i++) {
      cell_id[0] = i;
      grid->InsertNextCell(VTK_VERTEX, 1, cell_id);
    }
    idd->SetGrid(grid);
  } else {
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::SafeDownCast(idd->GetGrid());
    if(!grid) {
      vtkGenericWarningMacro("no adaptor grid to attach data to");
      return;
    }
    if(grid->GetNumberOfPoints() != static_cast<vtkIdType>(*count)) {
      grid->GetPoints()->SetNumberOfPoints(static_cast<vtkIdType>(*count));
      grid->Allocate(grid->GetNumberOfPoints());
      vtkIdType cell_id[1];
      for(vtkIdType i = 0; i < static_cast<vtkIdType>(*count); i++) {
        cell_id[0] = i;
        grid->InsertNextCell(VTK_VERTEX, 1, cell_id);
      }
      for(int i = 0; i < grid->GetPointData()->GetNumberOfArrays(); i++) {
        grid->GetPointData()->GetArray(i)->SetNumberOfTuples(static_cast<vtkIdType>(*count));
      }
    }
  }
}

void set_particle_position(long long int* id, double* x, double* y, double* z, const char* grid_name) {
  if(!vtkCPPythonAdaptorAPI::GetCoProcessorData()) {
    vtkGenericWarningMacro("unable to access coprocessor data");
    return;
  }
  vtkSmartPointer<vtkCPInputDataDescription> idd = vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName(grid_name);
  if(!idd) {
    vtkGenericWarningMacro("cannot access data description");
    return;
  }
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::SafeDownCast(idd->GetGrid());
  if(!grid) {
    vtkGenericWarningMacro("no adaptor grid to attach data to");
    return;
  }
  grid->GetPoints()->SetPoint(static_cast<vtkIdType>(*id), *x, *y, *z);
}

void set_scalar_particle(long long int* id, double* scalar, char* name, const char* grid_name) {
  if(!vtkCPPythonAdaptorAPI::GetCoProcessorData()) {
    vtkGenericWarningMacro("unable to access coprocessor data");
    return;
  }
  vtkSmartPointer<vtkCPInputDataDescription> idd = vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName(grid_name);
  if(!idd) {
    vtkGenericWarningMacro("cannot access data description");
    return;
  }
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::SafeDownCast(idd->GetGrid());
  if(!grid) {
    vtkGenericWarningMacro("no adaptor grid to attach data to");
    return;
  }
  if(grid->GetPointData()->GetArray(name) == nullptr) {
    vtkDoubleArray* array = vtkDoubleArray::New();
    array->SetName(name);
    array->SetNumberOfComponents(1);
    array->SetNumberOfTuples(grid->GetNumberOfPoints());
    //array->SetVoidArray(px.data(), static_cast<vtkIdType>((*components) * grid->GetNumberOfPoints()), 1);
    grid->GetPointData()->AddArray(array);
    array->Delete();
  }
  grid->GetPointData()->GetArray(name)->SetTuple1(static_cast<vtkIdType>(*id), *scalar);
}

void set_vector_particle(long long int* id, double* x_comp, double* y_comp, double* z_comp, char* name, const char* grid_name) {
  if(!vtkCPPythonAdaptorAPI::GetCoProcessorData()) {
    vtkGenericWarningMacro("unable to access coprocessor data");
    return;
  }
  vtkSmartPointer<vtkCPInputDataDescription> idd = vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName(grid_name);
  if(!idd) {
    vtkGenericWarningMacro("cannot access data description");
    return;
  }
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::SafeDownCast(idd->GetGrid());
  if(!grid) {
    vtkGenericWarningMacro("no adaptor grid to attach data to");
    return;
  }
  if(grid->GetPointData()->GetArray(name) == nullptr) {
    vtkDoubleArray* array = vtkDoubleArray::New();
    array->SetName(name);
    array->SetNumberOfComponents(3);
    array->SetNumberOfTuples(grid->GetNumberOfPoints());
    //array->SetVoidArray(px.data(), static_cast<vtkIdType>((*components) * grid->GetNumberOfPoints()), 1);
    grid->GetPointData()->AddArray(array);
    array->Delete();
  }
  grid->GetPointData()->GetArray(name)->SetTuple3(static_cast<vtkIdType>(*id), *x_comp, *y_comp, *z_comp);
}

void set_scalar_field(int* rank, double* scalar, char* name, const char* grid_name) {
  if(!vtkCPPythonAdaptorAPI::GetCoProcessorData()) {
    vtkGenericWarningMacro("unable to access coprocessor data");
    return;
  }
  vtkSmartPointer<vtkCPInputDataDescription> idd = vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName(grid_name); 
  if(!idd) {
    vtkGenericWarningMacro("cannot access data description");
    return;
  }
  vtkSmartPointer<vtkImageData> grid = vtkImageData::SafeDownCast(idd->GetGrid());
//  vtkSmartPointer<vtkMultiBlockDataSet> multi_block = vtkMultiBlockDataSet::SafeDownCast(idd->GetGrid());
//  vtkSmartPointer<vtkMultiPieceDataSet> multi_piece = vtkMultiPieceDataSet::SafeDownCast(multi_block->GetBlock(0));
//  vtkSmartPointer<vtkRectilinearGrid> grid = vtkRectilinearGrid::SafeDownCast(multi_piece->GetPiece(*rank));
  if(!grid) {
    vtkGenericWarningMacro("no adaptor grid to attach data to");
    return;
  }
  if(grid->GetPointData()->GetArray(name) == nullptr) {
    vtkSOADataArrayTemplate<double>* field = vtkSOADataArrayTemplate<double>::New();
    field->SetName(name);
    field->SetNumberOfComponents(1);
    field->SetNumberOfTuples(grid->GetNumberOfPoints());
    field->SetArray(0, scalar, grid->GetNumberOfPoints(), false, true);
    grid->GetPointData()->AddArray(field);
    field->Delete();
  }
  //if(update_field) {
  //  vtkSmartPointer<vtkDoubleArray> array = vtkDoubleArray::SafeDownCast(grid->GetPointData()->GetArray(name));
  //  array->SetArray(input_field, static_cast<vtkIdType>((*components) * grid->GetNumberOfPoints()), 1);
  //}
}

void set_vector_field(int* rank, double* x_comp, double* y_comp, double* z_comp, char* name, const char* grid_name) {
 if(!vtkCPPythonAdaptorAPI::GetCoProcessorData()) {
    vtkGenericWarningMacro("unable to access coprocessor data");
    return;
  }
  vtkSmartPointer<vtkCPInputDataDescription> idd = vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName(grid_name);
  if(!idd) {
    vtkGenericWarningMacro("cannot access data description");
    return;
  }
  vtkSmartPointer<vtkImageData> grid = vtkImageData::SafeDownCast(idd->GetGrid());
//  vtkSmartPointer<vtkMultiBlockDataSet> multi_block = vtkMultiBlockDataSet::SafeDownCast(idd->GetGrid());
//  vtkSmartPointer<vtkMultiPieceDataSet> multi_piece = vtkMultiPieceDataSet::SafeDownCast(multi_block->GetBlock(0));
//  vtkSmartPointer<vtkRectilinearGrid> grid = vtkRectilinearGrid::SafeDownCast(multi_piece->GetPiece(*rank));
  if(!grid) {
    vtkGenericWarningMacro("no adaptor grid to attach data to");
    return;
  }
  if(grid->GetPointData()->GetArray(name) == nullptr) {
    vtkSOADataArrayTemplate<double>* field = vtkSOADataArrayTemplate<double>::New();
    field->SetName(name);
    field->SetNumberOfComponents(3);
    field->SetNumberOfTuples(grid->GetNumberOfPoints());
    field->SetArray(0, x_comp, grid->GetNumberOfPoints(), false, true);
    field->SetArray(1, y_comp, grid->GetNumberOfPoints(), false, true);
    field->SetArray(2, z_comp, grid->GetNumberOfPoints(), false, true);
    grid->GetPointData()->AddArray(field);
    field->Delete();
  }
  //if(update_field) {
  //  vtkSmartPointer<vtkDoubleArray> array = vtkDoubleArray::SafeDownCast(grid->GetPointData()->GetArray(name));
  //  array->SetArray(input_field, static_cast<vtkIdType>((*components) * grid->GetNumberOfPoints()), 1);
  //}
}

add_executable (sz_test sz_compress_test.cpp)
target_include_directories (sz_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries (sz_test cpsz)

add_executable (sz_test_with_eb sz_compress_test_with_eb.cpp)
target_include_directories (sz_test_with_eb PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries (sz_test_with_eb cpsz)

add_executable (sz_compress_cp_preserve_2d_test sz_compress_cp_preserve_2d_test.cpp)
target_include_directories (sz_compress_cp_preserve_2d_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries (sz_compress_cp_preserve_2d_test cpsz)

add_executable (sz_compress_cp_preserve_3d_test sz_compress_cp_preserve_3d_test.cpp)
target_include_directories (sz_compress_cp_preserve_3d_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries (sz_compress_cp_preserve_3d_test cpsz FTK::libftk)

add_executable (sz_compress_cp_preserve_3d_unstructured_test sz_compress_cp_preserve_3d_unstructured_test.cpp)
target_include_directories (sz_compress_cp_preserve_3d_unstructured_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries (sz_compress_cp_preserve_3d_unstructured_test cpsz)

add_executable (sz_compress_2d_with_eb sz_compress_2d_with_eb.cpp)
target_include_directories (sz_compress_2d_with_eb PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries (sz_compress_2d_with_eb cpsz)

add_executable (cp_extraction_2d cp_extraction_2d.cpp)
target_include_directories (cp_extraction_2d PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries (cp_extraction_2d cpsz FTK::libftk)

add_executable (cp_extraction_3d cp_extraction_3d.cpp)
target_include_directories (cp_extraction_3d PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries (cp_extraction_3d cpsz FTK::libftk)

install (TARGETS sz_compress_cp_preserve_2d_test sz_compress_cp_preserve_3d_test sz_compress_cp_preserve_3d_unstructured_test sz_compress_2d_with_eb
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)
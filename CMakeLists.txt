if(DEFINED ENV{IDF_PATH})
    file(GLOB srcs
        "src/*.c"
        "src/*.cpp"
    )

    idf_component_register(
        SRCS ${srcs}
        INCLUDE_DIRS "./include"
    )
else()
    # Linux环境下的配置
    cmake_minimum_required(VERSION 3.10)
    project(rs_fec)

    # 设置C标准
    set(CMAKE_C_STANDARD 99)
    set(CMAKE_C_STANDARD_REQUIRED ON)

    # 源文件
    set(SRC_FILES
        "src/fec.cpp"
    )

    # 头文件目录
    include_directories(include)

    # 生成静态库
    add_library(rs_fec STATIC ${SRC_FILES})

    # 安装配置
    install(TARGETS rs_fec
        ARCHIVE DESTINATION lib
    )
    
    install(DIRECTORY include/
        DESTINATION include/rs_fec
    )
endif()
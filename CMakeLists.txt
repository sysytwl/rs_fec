file(GLOB srcs
    "src/*.c"
    "src/*.cpp"
)

idf_component_register(
    SRCS ${srcs}
    INCLUDE_DIRS "./include"
)

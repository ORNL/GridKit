# -*- Python -*-

with section("format"):

    # How wide to allow formatted cmake files
    line_width = 100

    # How many spaces to tab for indent
    tab_size = 2

    # If true, separate flow control names from their parentheses with a space
    separate_ctrl_name_with_space = False

    # If true, separate function names from parentheses with a space
    separate_fn_name_with_space = False

    # If a statement is wrapped to more than one line, then dangle the closing
    # parenthesis on its own line.
    dangle_parens = False

    # If a positional argument group contains more than this many arguments, then
    # force it to a vertical layout.
    max_pargs_hwrap = 3

    # Force vertical alignment if command takes more than this number of lines
    max_rows_cmdline = 1

    # List of command names which should always be wrapped
    always_wrap = ['configure_package_config_file',
                   'list',
                   'set_target_properties',
                   'target_include_directories', 
                   'target_link_libraries',
                   'FILES',
                   'HEADERS',
                   'INCLUDE_DIRECTORIES',
                   'LINK_LIBRARIES',
                   'NAMES',
                   'PATHS',
                   'PATH_SUFFIXES',
                   'SOURCES'
                   ]
    
with section("parse"):

    # Formatting for custom macros
    additional_commands = {
        "gridkit_add_library": {
            "pargs": 1,  # Number of initial positional arguments
            "flags": [],
            "kwargs": {
                "SOURCES": "*",
                "HEADERS": "*",
                "LINK_LIBRARIES": "*",
                "INCLUDE_DIRECTORIES": "*",
                "COMPILE_OPTIONS": "*",
            }
        },
        "enzyme_build_object": {
            "flags": [],
            "kwargs": {
                "NAME": 1,
                "SOURCES": "*",
                "LINK_LIBRARIES": "*",
                "INCLUDE_DIRECTORIES": "*",
            }
        }
    }

with section("markup"):

    # Globally disable comment markup processing
    enable_markup = False

    # Do not reflow the first comment block in each file
    first_comment_is_literal = True

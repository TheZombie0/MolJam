import sys


def main():
    if "--mcp" in sys.argv:
        from .mcp_server import main as mcp_main

        mcp_main()
    else:
        print("moljam - Molecular database quality scoring toolkit")
        print("  python -m moljam --mcp    Start MCP server")


if __name__ == "__main__":
    main()

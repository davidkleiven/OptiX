import sys
sys.path.append("Applications/WgBuilder")
import editor as ed
import tkinter as tk


def main():
    root = tk.Tk()
    root.wm_title("Wave Guide Builder")
    editor = ed.Editor( root )
    root.mainloop()

if __name__ == "__main__":
    main()

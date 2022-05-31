import tkinter as tk
from tkinter import filedialog

from fonction import function


root = tk.Tk()
root.title('Qualité NGS')
root.iconbitmap('logo-polytech.ico')
root.geometry("1000x800")

apps = []
results = []
def addApp():
    filename = filedialog.askopenfilename(initialdir="/", title="Select File", multiple=True)
    apps.extend(filename)
    for widget in frame.winfo_children():
        widget.destroy()   
    
    for app in apps:
        label = tk.Label(frame, text=app, bg="white")
        label.pack()

def runApps():
   
    res_text = ["Classe de la qualité", "Probabilités d'appartenance pour chaque classe", "Pourcentage des séquences OPT"]
    df_amp = function(apps)
    if len(df_amp) !=0 :    
        for i in range(len(res_text)):
            labelF = tk.LabelFrame(frame, text=res_text[i], bg="white")
            labelF.pack(padx=10, pady=10)
            
            label = tk.Label(labelF, text = df_amp[i], bg="white", )
            label.pack()
    else :
        label = tk.Label(frame, text = "Qualité NC. Le fichier Variant Caller n'a aucune mutation ayant dépassé le seuil de qualité.", bg="white")
        label.pack(padx=10, pady=10) 
      
        
canvas = tk.Canvas(root, bg="#263D42")
canvas.place(relwidth=0.8, relheight=0.8, relx=0.1, rely=0.1)

frame = tk.Frame(canvas, bg="white")
frame.place(relwidth=0.8, relheight=0.8, relx=0.1, rely=0.1)

openApps = tk.Button(root, text="Run", padx=10, pady=5, fg="white", bg="#263D42", command=runApps, cursor="hand2")
openApps.pack(side="bottom")

openFile = tk.Button(root, text="Open File", padx=10, pady=5, fg="white", bg="#263D42", command=addApp, cursor="hand2")
openFile.pack(side="bottom")



root.mainloop()
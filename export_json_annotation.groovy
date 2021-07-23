import com.google.gson.Gson
import qupath.lib.io.GsonTools
import qupath.lib.objects.PathObject


import qupath.lib.objects.PathObject

def imageData = QPEx.getCurrentImageData()
def server = imageData.getServer()
String FilePath = server.getFile()
print FilePath

String[] parts = FilePath.split("/")
String imageName = parts[-1]
print imageName
String[] parts2 = imageName.split("\\.")
String imageNumber = parts2[0]
print imageNumber
//================================================================
//OutputPath is the json file you want to save to. Please use'/' instead of'\'!!!
//------------Please enter your own json OutputPath:-----------------
String OutputPath="/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/"
String suffix=".json"
String pathCompleate=OutputPath+imageNumber+suffix
print pathCompleate

def path = buildFilePath(pathCompleate)
def file = new File(path)
file.text = ''
def annotations = getAnnotationObjects()
//boolean prettyPrint = true
def gson = GsonTools.getInstance(false)
def text=gson.toJson(annotations)
println text
//file << text

file.setText(text, 'UTF-8')
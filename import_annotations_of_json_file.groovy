import com.google.gson.Gson;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonParser;

public static String readFile(String filename) {
        String result = "";
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            StringBuilder sb = new StringBuilder();
            String line = br.readLine();
            while (line != null) {
                sb.append(line);
                line = br.readLine();
            }
            result = sb.toString();
        } catch(Exception e) {
            e.printStackTrace();
        }
        return result;
    }

import com.google.gson.Gson
import qupath.lib.io.GsonTools
import qupath.lib.objects.PathObject

def imageData = QPEx.getCurrentImageData()
def server = imageData.getServer()
String FilePath = server.getPath()
print FilePath

String[] parts = FilePath.split("/")
String imageName = parts[-1]
print imageName
String[] parts2 = imageName.split("\\.")
String imageNumber = parts2[0]
print imageNumber
//================================================================
//InputPath is the json file belongs to. Please use'/' instead of'\'!!!
//Such as : "C:/Users/10158/Desktop/XI/"
//------------Please enter your own json InputPath:-----------------
String InputPath='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/05_beijing/02_qupath_res/'
//"/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/02_qupath_res/"
//"D:/00_Apr05_CRC_Liver_Image/06_Corrected_image_json/"

String suffix=".json"
String pathCompleate=InputPath+imageNumber+suffix
print pathCompleate
JsonElement je = new JsonParser().parse(readFile(pathCompleate))
print(je)


import com.google.gson.Gson
import qupath.lib.io.GsonTools
import qupath.lib.objects.PathAnnotationObject
boolean prettyPrint = true
def gson = GsonTools.getInstance(prettyPrint)

JsonArray ja = je.getAsJsonArray();
Iterator itr = ja.iterator();
            while (itr.hasNext()) {
                JsonElement je1 = (JsonElement) itr.next();
                def pathObject2 = gson.fromJson(je1, PathAnnotationObject)
                print(pathObject2)
                addObject(pathObject2)
            }


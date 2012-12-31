function favBrowser() {
    var list=document.getElementById("bandList");
    document.getElementById("choice").value=list.options[list.selectedIndex].text;
    }
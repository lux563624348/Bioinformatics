
var cfo = {
	getDataFromMetaTags: function() {
		var data = {   vin: this.fromMetaTag('cfo:vin'),
					   accountname: this.fromMetaTag('cfo:accountname'),
					   password: this.fromMetaTag('cfo:password'),
					   eventType: this.fromMetaTag('cfo:eventType'),
					   suggest: this.fromMetaTag('cfo:suggest'),
					   vinSuggest: this.fromMetaTag('cfo:vinSuggest')
		};

		if(data.suggest == '') {
		   data.suggest = data.vinSuggest;
		}
		return data;
	},

	upgradeVhrToCip: function() {
		var cipData = this.getDataFromMetaTags();

		if(cipData.eventType.indexOf("PHLEligible") == -1) {
		   cipData.eventType = "CIPUpgradeFromVHR";
		}
		else {
		   cipData.eventType = "CIPUpgradeFromVHR-PHLEligible";
		}

		doPost('/cfm/cfoEventHandler.cfm', cipData);
	},

 	toggleToLanguage: function(language) {
		// Language values should be ES or EN
		if(language != "es" && language != "en")
			return;

		var toggleData = this.getDataFromMetaTags();
		toggleData.language = language;
		toggleData.skipSignature='Y';
		toggleData.eventType = "RunReport" + this.fromMetaTag('cfo:reportType');

		doPost('/cfm/cfoEventHandler.cfm', toggleData);
	},

    fromMetaTag: function(name){
		var value = $('meta[name="'+name+'"]').attr('content');
		return (value) ? value : "";
	}
};

function doPost(url, data) {
	var form = "<form action='"+url+"' method='post' name='generatedForm'>";
	for(var key in data){
		form += "<input type='hidden' name='"+key+"' value='"+data[key]+"'/>"
	}
	form +="</form>";

	var innerHtmlMaker = document.createElement('span');
	innerHtmlMaker.innerHTML = form;

	document.body.appendChild(innerHtmlMaker.firstChild);

	document.forms['generatedForm'].submit();
}